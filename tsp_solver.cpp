#include "tsp_solver.h"

#define PSEP_PIV_FRAC 0
#define PSEP_PIV_SUBTOUR 1
#define PSEP_PIV_TOUR 2
#define PSEP_PIV_OPT_TOUR 3

using namespace std;

bool TSP_Solver::devex = true;

int TSP_Solver::pivot_until_change(int *pivot_status_p){
  int rval = 0;
  int itcount = 0, icount = 0;
  int rowcount = PSEPlp_numrows(&m_lp);
  bool integral = false, conn = false, dual_feas = false;
  *pivot_status_p = -1;

  double routine_start, round_start = PSEP_zeit();

  old_rowstat.resize(rowcount);

  while(true){
    if(++itcount == 3 * m_graph.node_count){
      rval = 1;
      cerr << "Pivot terminated due to iteration limit" << endl;
      goto CLEANUP;
    }

    if((dual_feas = is_dual_feas()))
      break;

    routine_start = PSEP_zeit();
    rval = PSEPlp_getbase(&m_lp, &old_colstat[0], &old_rowstat[0]);
    if(rval) goto CLEANUP;
    PSEP_Timer::getbase_time += PSEP_zeit() - routine_start;

    rval = primal_pivot();
    if(rval) goto CLEANUP;

    rval = set_edges();
    if(rval) goto CLEANUP;

    if(fabs(get_obj_val() - m_min_tour_value) >= LP_EPSILON)
      break;    
  }

  round_start = PSEP_zeit() - round_start;

  cout << itcount << " pivots performed in "
       << round_start << "s, pivot obj val: " << get_obj_val()
       << endl;

  if(!devex && (itcount > m_graph.node_count / 2)){
    rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_Simplex_PGradient,
			CPX_PPRIIND_DEVEX);
    if(rval)
      goto CLEANUP;
    else{
      cout << "////////Switched to devex/////////" << endl;
      devex = true;
    }
  }

  rval = set_support_graph();
  if(rval) goto CLEANUP;

  integral = is_integral();
    if(integral){
    conn = G_Utils::connected(&G_s, &icount, island, 0);
    if(integral && conn){
      if(dual_feas)
	*pivot_status_p = PSEP_PIV_OPT_TOUR;
      else
	*pivot_status_p = PSEP_PIV_TOUR;
    } else
      *pivot_status_p = PSEP_PIV_SUBTOUR;
  } else
    *pivot_status_p = PSEP_PIV_FRAC;

 CLEANUP:
    if(rval)
      cerr << "Error entry point: pivot_until_change" << endl;
    return rval;
}

int TSP_Solver::pivot_back(){
  cout << "Calling pivot back...";

  int rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);

  if(rval)
    cerr << "Error entry point: pivot_back()" << endl;
  else
    cout << "Pivoted back successfully" << endl;
  return rval;
}


//Return false if not a tour
bool TSP_Solver::update_current_tour_indices(vector<int> &tour_nodes){
  int icount = 0;
  if(!G_Utils::connected(&G_s, &icount, island, 0)){
    cout << "Best tour is not a connected graph" << endl;
    return false;
  }

  for(int i = 0; i < m_graph.node_count; i++){
    tour_nodes[i] = island[i];
  }

  return true;
}

int TSP_Solver::update_best_tour(){
  int rval = 0;
  double objval = 0;

  rval = (!update_current_tour_indices(best_tour_nodes));
  
  if (rval){
    fprintf(stderr, "New best tour is not a tour!!\n");
    goto CLEANUP;
  }

  for (int i = 0; i < m_graph.edge_count; i++){
    if (m_lp_edges[i] < LP_EPSILON){
      best_tour_edges[i] = 0;
    } else {
      best_tour_edges[i] = 1;
      objval += m_graph.edges[i].len;
    }
  }

  if(objval > m_min_tour_value){
    fprintf(stderr, "New best tour is worse!\n");
    rval = 1;
    goto CLEANUP;
  }

  m_min_tour_value = objval;

  for(int i = 0; i < m_graph.node_count; i++){
    perm[best_tour_nodes[i]] = i;
  }

 CLEANUP:
  return rval;
}
