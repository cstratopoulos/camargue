#include "tsp_solver.h"

using namespace std;

int TSP_Solver::basis_init(){  
  int rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);
  if(rval) goto CLEANUP;

  rval = primal_pivot();
  if(rval) goto CLEANUP;
    
  rval = set_edges();
  if(rval) goto CLEANUP;

  if(fabs(get_obj_val() - m_min_tour_value) >= LP::EPSILON){
    cerr << "BASIS INIT SWITCHED OBJ VAL!!!! TO: "
	 << get_obj_val() << endl;
    return 1;
  }
  
  for(int i = 0; i < m_graph.edge_count; i++){
    if(fabs(m_lp_edges[i] - best_tour_edges[i]) >= LP::EPSILON){
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


bool TSP_Solver::is_dual_feas(){
  return PSEPlp_dualfeas(&m_lp);
}

bool TSP_Solver::is_integral(){
  bool result = true;
  for(int i = 0; i < m_lp_edges.size(); i++){
    if ((m_lp_edges[i] >= LP::EPSILON) && (m_lp_edges[i] <= (1 - LP::EPSILON))){
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
    if(m_lp_edges[i] >= LP::EPSILON){
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
