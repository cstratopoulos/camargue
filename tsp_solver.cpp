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
  int rowcount = PSEPlp_numrows(&m_lp);
  bool found_edge_change = false, found_bhead_change = false,
    integral = false, conn = false, dual_feas = false;
  *old_nb_stat_p = -1;
  *pivot_status_p = -1;

  basis_headers[0].resize(rowcount);
  basis_headers[1].resize(rowcount);

  while(true){
    if(++itcount == 3 * m_graph.node_count){
      rval = 1;
      cerr << "Pivot terminated due to iteration limit" << endl;
      goto CLEANUP;
    }

    if((dual_feas = is_dual_feas()))
      break;

    rval = (PSEPlp_bhead(&m_lp, &old_header[0], NULL) ||
	    PSEPlp_getbase(&m_lp, &old_base[0], NULL));
    if(rval) goto CLEANUP;
      

    rval = primal_pivot();
    if(rval) goto CLEANUP;

    rval = set_edges();
    if(rval) goto CLEANUP;

    rval = PSEPlp_bhead(&m_lp, &new_header[0], NULL);
    if(rval) goto CLEANUP;

    if(fabs(get_obj_val() - m_min_tour_value) >= LP_EPSILON)
      break;

    for(int i = 0; i < rowcount; i++){
      if(old_header[i] != new_header[i]){
	found_bhead_change = true;

	if(new_header[i] >= 0){//a nonbasic edge was pivoted into basis
	  if(fabs(m_lp_edges[new_header[i]] - best_tour_edges[new_header[i]])
	     >= LP_EPSILON){
	    found_edge_change = true;

	    *old_b_p = old_header[i];
	    *old_nb_p = new_header[i];
	    *old_nb_stat_p = old_base[new_header[i]];

	    break;
	  }
	}

	if(old_header[i] >= 0){//newheader[i] < 0, i.e., a slack was pivoted in
	  if(fabs(m_lp_edges[old_header[i]] - best_tour_edges[new_header[i]])
	     >= LP_EPSILON){
	    found_edge_change = true;

	    *old_b_p = old_header[i];
	    *old_nb_p = new_header[i];

	    int row_index = - new_header[i] - 1;
	    char sense;

	    rval = PSEPlp_getsense(&m_lp, &sense, row_index);
	    if(rval) goto CLEANUP;

	    switch(sense){
	    case 'G':
	      *old_nb_stat_p = CPX_AT_LOWER;
	      break;
	    case 'L':
	      *old_nb_stat_p = CPX_AT_UPPER;
	      break;
	    default:
	      cerr << "Problem getting sense of slack row" << endl;
	      rval = 1;
	      goto CLEANUP;
	    }

	    break;
	  }
	}
      }
    }

    if(found_edge_change)
      break;

    if(!found_bhead_change){
      for(int i = 0; i < ecount; i++)
	if(fabs(m_lp_edges[i] - best_tour_edges[i]) >= LP_EPSILON){
	  *old_b_p = i;
	  *old_nb_p = i;
	  *old_nb_stat_p = old_base[i];
	  
	  found_edge_change = true;
	  break;
	}
    }

    if(found_edge_change)
      break;
  }

  cout << itcount << " pivots performed, pivot obj val: " << get_obj_val()
       << endl;

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
