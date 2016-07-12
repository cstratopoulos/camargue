#include "LPcore.h"

#include <math.h>
#include <iomanip>

using namespace std;

inline bool PSEP_LP_Core::is_dual_feas(){
  return PSEPlp_dualfeas(&m_lp);
}

inline bool PSEP_LP_Core::is_integral(){
  for(int i = 0; i < m_lp_edges.size(); i++)
    if((m_lp_edges[i] >= LP::EPSILON) && (m_lp_edges[i] <= (1 - LP::EPSILON)))
      return false;

  return true;
}

int PSEP_LP_Core::pivot(){
  int infeasible = 0;
  int rval = PSEPlp_primal_pivot(&m_lp, &infeasible);

  if(rval)
    cerr << "Entry point LP_Core::pivot(), infeasible " << infeasible << "\n";
  return rval;
}


int PSEP_LP_Core::pivot_back(){
  int rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);

  if(rval)
    cerr << "Error entry point: pivot_back()" << endl;
  return rval;
}

double PSEP_LP_Core::get_obj_val(){
  double objval;
  PSEPlp_objval(&m_lp, &objval);
  return objval;
}

double PSEP_LP_Core::set_edges(){
  int rval = 0;

  if (m_lp_edges.size() < m_graph.edge_count)
    m_lp_edges.resize(m_graph.edge_count);

  rval = PSEPlp_x (&m_lp, &m_lp_edges[0]);
  if(rval)
    fprintf(stderr, "failed to set_edges(), rval %d\n", rval);

  return rval;
}

int PSEP_LP_Core::basis_init(){
  int rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);
  if(rval) goto CLEANUP;

  rval = pivot();
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

double PSEP_LP_Core::set_support_graph(){
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

int PSEP_LP_Core::update_best_tour(){
  double objval = 0;
  
  for(int i = 0; i < m_graph.node_count; i++)
    best_tour_nodes[i] = island[i];

  for(int i = 0; i < m_graph.edge_count; i++)
    if(m_lp_edges[i] < LP::EPSILON)
      best_tour_edges[i] = 0;
    else {
      best_tour_edges[i] = 1;
      objval += m_graph.edges[i].len;
    }

  if(objval > m_min_tour_value){
    cerr << "New best tour is worse!\n";
    return 1;
  }

  m_min_tour_value = objval;

  for(int i = 0; i < m_graph.node_count; i++)
    perm[best_tour_nodes[i]] = i;

  return 0;
};

int PSEP_LP_Core::pivot_until_change(int *pivot_status_p){
int rval = 0;
  int itcount = 0, icount = 0;
  int rowcount = PSEPlp_numrows(&m_lp);
  bool integral = false, conn = false, dual_feas = false;
  *pivot_status_p = -1;

  double round_start = PSEP_zeit();

  old_rowstat.resize(rowcount);

  while(true){
    if(++itcount == 3 * rowcount)
      if(LP::PRICING::BEHAVIOR::choice = LP::PRICING::BEHAVIOR::DYNAMIC)
	change_pricing();

    if((dual_feas = is_dual_feas()))
      break;

    rval = PSEPlp_getbase(&m_lp, &old_colstat[0], &old_rowstat[0]);
    if(rval) goto CLEANUP;

    rval = pivot();
    if(rval) goto CLEANUP;

    rval = set_edges();
    if(rval) goto CLEANUP;

    if(fabs(get_obj_val() - m_min_tour_value) >= LP::EPSILON)
      break;    
  }

  round_start = PSEP_zeit() - round_start;

  rval = set_support_graph();
  if(rval) goto CLEANUP;

  integral = is_integral();
    if(integral){
    conn = G_Utils::connected(&G_s, &icount, island, 0);
    if(integral && conn){
      if(dual_feas)
	*pivot_status_p = PIVOT::FATHOMED_TOUR;
      else
	*pivot_status_p = PIVOT::TOUR;
    } else
      *pivot_status_p = PIVOT::SUBTOUR;
  } else
      *pivot_status_p = PIVOT::FRAC;

    cout << "Did " << itcount << " pivots in "
	 << setprecision(2) << round_start << "s, "
	 << "obj val: " << setprecision(6) << get_obj_val() << "\n";

 CLEANUP:
    if(rval)
      cerr << "Error entry point: pivot_until_change" << endl;
    return rval;
}

void PSEP_LP_Core::change_pricing(){
  int newprice;

  cout << "//////SWITCHED PRICING TO ";
  switch(LP::PRICING::choice){
  case LP::PRICING::DEVEX:
    newprice = CPX_PPRIIND_DEVEX;
    cout << "DEVEX/////\n";
    break;
  case LP::PRICING::STEEPEST:
    newprice = CPX_PPRIIND_STEEPQSTART;
    cout << "STEEPEST EDGE W SLACK INITIAL NORMS\n";
    break;
  case LP::PRICING::STEEPEST_REAL:
    newprice = CPX_PPRIIND_STEEP;
    cout << "GENUINE STEEPEST EDGE////\n";
    break;
  }

  if(CPXsetintparam(lp->cplex_env, CPXPARAM_Simplex_PGradient, newprice))
    cerr << "ERROR: PRICING SWITCH DID NOT TAKE PLACE\n";
  LP::PRICING::SWITCHING::choice = LP::PRICING::BEHAVIOR::OFF;
}

