#include "tsp_solver.h"

#define PSEP_PIV_FRAC 0
#define PSEP_PIV_SUBTOUR 1
#define PSEP_PIV_TOUR 2
#define PSEP_PIV_OPT_TOUR 3

using namespace std;

int TSP_Solver::pivot_until_change(int *old_b_p, int *old_nb_p,
				   int *old_nb_stat_p, int *pivot_status_p){
  int rval = 0;
  int ecount = m_graph.edge_count, itcount = 0, icount = 0;
  bool found_change = false;
  bool integral = false, conn = false;
  *old_b_p = -1;
  *old_nb_p = -1;
  *old_nb_stat_p = -1;
  *pivot_status_p = PSEP_PIV_TOUR;

  int rowcount = PSEPlp_numrows(&m_lp);
  vector<int> oldhead(rowcount); vector<int> newhead(rowcount);

  PSEPlp_bhead(&m_lp, &oldhead[0], NULL);

  vector<int>rowbase(rowcount);

  bool dual_feas = false;

  while(true){
    if((rval = (++itcount == 3 * m_graph.node_count))){
      cerr << "Pivot terminated due to iteration limit" << endl;
      goto CLEANUP;
    }

    dual_feas = is_dual_feas();
    if(dual_feas){
      //cout << "Pivot is optimal" << endl;
      break;
    }

    rval = primal_pivot();
    if(rval) goto CLEANUP;




    rval = set_edges();
    if(rval) goto CLEANUP;

    
    if(fabs(get_obj_val() - m_min_tour_value) >= LP_EPSILON)
      found_change = true;    
    else{
      for(int i = 0; i < ecount; i++){
	if(fabs(m_lp_edges[i] - best_tour_edges[i]) < LP_EPSILON)
	  continue;

	found_change = true;
	break;
      }
    }
    

    if(found_change){
      PSEPlp_getbase(&m_lp, &new_base[0], NULL);
      PSEPlp_bhead(&m_lp, &newhead[0], NULL);
      break;
      
    }

    PSEPlp_getbase(&m_lp, &old_base[0], NULL);
    PSEPlp_bhead(&m_lp, &oldhead[0], NULL);
  }


  cout << itcount << " pivots performed, "
       << "Pivot obj val: " << get_obj_val() << endl;

  PSEPlp_getbase(&m_lp, NULL, &rowbase[0]);
    cout << "Checking conjecture on rowbase...";
    for(int i = m_graph.node_count; i < rowcount; i++){
      cout << "row " << i << " has status: "
	   << rowbase[i] << endl;
    }
    cout << "Done." << endl;

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
      goto CLEANUP;
    }

    *pivot_status_p = PSEP_PIV_SUBTOUR;
  } else
    *pivot_status_p = PSEP_PIV_FRAC;


  for(int i = 0; i < ecount; i++){
    if(old_base[i] != new_base[i] &&
       fabs(m_lp_edges[i] - best_tour_edges[i]) >= LP_EPSILON){
      if(old_base[i] + new_base[i] == 2){ //variable was moved from lower->upper
	//(0 -> 2)
	*old_b_p = i;
	*old_nb_p = i;
	*old_nb_stat_p = old_base[i];
	///*
	cout << "Basis change on column " << i
	     << ", status: " << old_base[i] << " -> " << new_base[i]
	     << ", solution vec: " << best_tour_edges[i] << " -> "
	     << m_lp_edges[i] << endl;//*/
	break;
      }

      if(old_base[i] == 1){
	*old_b_p = i;
	///*
	cout << "Basis change on column " << i
	     << ", status: " << old_base[i] << " -> " << new_base[i]
	     << ", solution vec: " << best_tour_edges[i] << " -> "
	     << m_lp_edges[i] << endl;
	//*/
      } else {
	*old_nb_p = i;
	*old_nb_stat_p = old_base[i];
	///*
	cout << "Basis change on column " << i
	     << ", status: " << old_base[i] << " -> " << new_base[i]
	     << ", solution vec: " << best_tour_edges[i] << " -> "
	     << m_lp_edges[i] << endl;
	//*/
      }

      if(*old_b_p != -1 && *old_nb_p != -1 && *old_nb_stat_p != -1)
	break;
    }
  }

  cout << "Now comparing oldhead and newhead for change..." << endl;
  for(int i = 0; i < rowcount; i++){
    if(oldhead[i] != newhead[i]){
      cout << "Disagreement: oldhead[" << i << "] = " << oldhead[i]
	   << ", newhead[" << i << "] = " << newhead[i] << endl;
    }
  }

  /*
  cout << "Now scanning for complete list of basis changes...." << endl;
  for(int i = 0; i < ecount; i++){
    if(old_base[i] != new_base[i]){
	cout << "Basis change on column " << i
	     << ", status: " << old_base[i] << " -> " << new_base[i]
	     << ", solution vec: " << best_tour_edges[i] << " -> "
	     << m_lp_edges[i] << endl;
    }
  }*/

  

 CLEANUP:
  if(rval)
    cerr << "Entry point: pivot_until_change" << endl;
  PSEPlp_getbase(&m_lp, &old_base[0], NULL);
  return rval;
}

int TSP_Solver::pivot_back(const int old_basic, const int old_nonbasic,
			   const int old_nb_stat){  
  cout << "Calling pivot back...";
  
  int rval = PSEPlp_pivot(&m_lp, old_basic, old_nonbasic, old_nb_stat);
  if(rval) goto CLEANUP;

  rval = set_edges();
  if(rval) goto CLEANUP;

  if(fabs(get_obj_val() - m_min_tour_value) >= LP_EPSILON){
    cerr << "Obj val disagreement, did not pivot back to tour!" << endl;
    cerr << "Obj val: " << get_obj_val() << ", min tour val: "
	 << m_min_tour_value << endl;
    rval = 1; goto CLEANUP;
  }

  for(int i = 0; i < m_graph.edge_count; i++){
    if(fabs(m_lp_edges[i] - best_tour_edges[i]) >= LP_EPSILON){
      cerr << "Entry disagreement, did not pivot back to tour!" << endl;
      rval = 1; goto CLEANUP;
    }
  }

 CLEANUP:
  if(rval)
    cerr << "Entry point: pivot_back" << endl;
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
