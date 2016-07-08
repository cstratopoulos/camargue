#include "tsp_solver.h"

using namespace std;

int TSP_Solver::basis_init(){
  vector<int> rowstat(m_graph.node_count);
  
  int rval = PSEPlp_copybase(&m_lp, &old_base[0], &rowstat[0]);
  if(rval) goto CLEANUP;

    
  rval = no_opt();
  if(rval) goto CLEANUP;

  rval = PSEPlp_bhead(&m_lp, &old_header[0], NULL);

  cout << "Comparing entries of basis header with basic status of columns:"
       << endl;
  for(int i = 0; i < basis_headers[0].size(); i++){
    cout << "Basis header " << i << ": " << old_header[i]
	 << ", colstat[" << old_header[i] << "] = "
	 << old_base[old_header[i]] << endl;
  }

  rval = primal_pivot();
  if(rval) goto CLEANUP;
    
  rval = set_edges();
  if(rval) goto CLEANUP;

  if(fabs(get_obj_val() - m_min_tour_value) >= LP_EPSILON){
    cerr << "BASIS INIT SWITCHED OBJ VAL!!!! TO: "
	 << get_obj_val() << endl;
    return 1;
  }
  
  for(int i = 0; i < m_graph.edge_count; i++){
    if(fabs(m_lp_edges[i] - best_tour_edges[i]) >= LP_EPSILON){
      cerr << "Found disagreement!!!!!: " << m_lp_edges[i] << " vs "
	   << best_tour_edges[i] << endl;
      return 1;
    }
  }
  
  
  CLEANUP:
  if(rval)
    cerr << "Error entry point: basis_init()" << endl;
  return rval;
}

void TSP_Solver::test_bhead(){
  cout << "Printing basic variables and corresponding LP solution values..."
       << endl;

  int numbasic = 0;
  for(int i = 0; i < m_graph.edge_count; i++){
    if(old_base[i] == CPX_BASIC){
      cout << "Variable " << i << " is basic, lp value: "
	   << m_lp_edges[i] << endl;
      numbasic++;
    }
  }

  cout << "Total number of variables in basis: " << numbasic << endl;

  int numrows = PSEPlp_numrows(&m_lp);
  vector<int> head(numrows);
  vector<double> x(numrows);

  cout << "Called bhead, rval: " << PSEPlp_bhead(&m_lp, &head[0], &x[0])
       << endl;

  cout << "Printing every entry of head and corresponding x..." << endl;

  for(int i = 0; i < numrows; i++){
    cout << "head[" << i << "] = " << head[i] << ", x[" << i << "] = "
	 << x[i];
    if(head[i] < 0){
      int row = (-head[i]) - 1;
      cout << ", lp_edges[" << row << "] = " << m_lp_edges[row]
	   << ", old_base[" << row << "] = " << old_base[row] << endl;
    } else {
      cout << ", lp_edges[" << head[i] << "] = " << m_lp_edges[head[i]]
	   << ", old_base[" << head[i] << "] = " << old_base[head[i]] << endl;
    }
  }

  cout << "Checking if any lp entry is somehow more than 1???" << endl;

  for(int i = 0; i < m_graph.edge_count; i++){
    if(m_lp_edges[i] > 1)
      cout << "m_lp_edges[" << i << "] = " << m_lp_edges[i] << endl;
  }

  cout << "Done. " << endl;
  
}


bool TSP_Solver::is_dual_feas(){
  return PSEPlp_dualfeas(&m_lp);
}

bool TSP_Solver::is_integral(){
  bool result = true;
  for(int i = 0; i < m_lp_edges.size(); i++){
    if ((m_lp_edges[i] >= LP_EPSILON) && (m_lp_edges[i] <= (1 - LP_EPSILON))){
      result = false;
      break;
    }
  }

  return result;
}


int TSP_Solver::primal_opt() {
	int infeasible = 0;
	int rval = PSEPlp_primal_opt (&m_lp, &infeasible);

    if (rval) {
        throw "primal_opt failed"; 
    }

    return infeasible;
}

int TSP_Solver::no_opt(){
  int rval = PSEPlp_no_opt (&m_lp);
  if(rval){
    throw "no_opt failed";
  }

  return rval;
}

int TSP_Solver::primal_pivot() {
  int infeasible = 0;
  double start = PSEP_zeit();
  int rval = PSEPlp_primal_pivot (&m_lp, &infeasible);

  if (rval)
    cerr << "Entry point primal_pivot(), infeasible "
	 << infeasible << endl;
  PSEP_Timer::pivot_time += PSEP_zeit() - start;
  return rval;
}

int TSP_Solver::dual_pivot() {
  int infeasible = 0;
  int rval = PSEPlp_dual_pivot (&m_lp, &infeasible);

  if (rval) {
    throw "dual_pivot failed";
  }

  return infeasible;
}

double TSP_Solver::get_obj_val() {
	double objval;
	PSEPlp_objval(&m_lp, &objval);
	return objval;
}


int TSP_Solver::set_edges(){
  int rval = 0;
  
  
  if (m_lp_edges.size() < m_graph.edge_count)
    m_lp_edges.resize(m_graph.edge_count);

  rval = PSEPlp_x (&m_lp, &m_lp_edges[0]);
  if(rval)
    fprintf(stderr, "failed to set_edges(), rval %d\n", rval);

  return rval;
}

int TSP_Solver::set_support_graph(){
  int rval = 0;
  
  support_indices.clear();
  support_elist.clear();
  support_ecap.clear();

  for(int i = 0; i < m_graph.edge_count; i++){
    if(m_lp_edges[i] >= LP_EPSILON){
      support_indices.push_back(i);
      support_ecap.push_back(m_lp_edges[i]);
      support_elist.push_back(m_graph.edges[i].end[0]);
      support_elist.push_back(m_graph.edges[i].end[1]);
    }
  }

  rval = G_Utils::build_s_graph(m_graph.node_count, support_indices.size(),
				m_graph.edges, support_indices, m_lp_edges,
				&G_s);
  
  if(rval)
    cerr << "Problem setting support graph" << endl;
  return rval;
}

void TSP_Solver::print_lp_edges(){
  for(int i = 0; i < m_lp_edges.size(); i++){
    if(m_lp_edges[i] >= LP_EPSILON){
    cout << "Edge number " << i << " "
	 <<"[" << m_graph.edges[i].end[0] << ", " << m_graph.edges[i].end[1]
	 << "] " << "lp vec val: " << m_lp_edges[i]
	 << " len: " << m_graph.edges[i].len << endl;
      }
  }
}

void TSP_Solver::print_best_tour_edges(){
  for(int i = 0; i < best_tour_edges.size(); i++){
    if (best_tour_edges[i] == 1){
      cout << "Edge number " << i << " "
	   <<"[" << m_graph.edges[i].end[0] << ", " << m_graph.edges[i].end[1]
	   << "]"
	   <<" len: " << m_graph.edges[i].len << endl;
    }
  }
}

void TSP_Solver::print_best_tour_nodes(){
  for(int i = 0; i < best_tour_nodes.size(); i++){
    cout << best_tour_nodes[i] << endl;
  }
}

void TSP_Solver::printrow(const int rownum, int *rmatind0, int *rmatind1){
  int rmatbeg[1], rmatind[2], rmatspace = 2, surplus = 0, nzcnt = 0;
  double rmatval[2];
  rmatbeg[0] = 0;

  PSEPlp_getrows(&m_lp, &nzcnt, rmatbeg, rmatind, rmatval, rmatspace, &surplus,
		 rownum, rownum);

  *rmatind0 = rmatind[0]; *rmatind1 = rmatind[1];

}
