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

int PSEP_LP_Core::primal_opt(){
  int infeasible = 0;
  double start = PSEP_zeit();
  int rval = PSEPlp_primal_opt(&m_lp, &infeasible);
  start = PSEP_zeit() - start;

  if(rval)
    cerr << "Entry point LP_Core::primal_opt(), infeasible "
	 << infeasible << "\n";
  else
    cout << "Primal optimized with obj val " << get_obj_val()
	 << " in " << start << "s.\n";

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

int PSEP_LP_Core::rebuild_basis(){
  int rval = 0, infeas = 0;
  int ecount = m_lp_edges.size();
  int ncount = best_tour_nodes.size();
  int num_removed = 0;
  double objval;
  vector<double> tour_obj(ecount);
  vector<double> old_obj(ecount);
  vector<int> lp_indices(ecount);
  old_colstat.resize(ecount);
  double prune;

  rval = PSEPlp_getobj(&m_lp, &old_obj[0], ecount);
  if(rval) return rval;

  for(int i = 0; i < ecount; i++){
    tour_obj[i] = 1 - (2 * best_tour_edges[i]);
    lp_indices[i] = i;
  }

  rval = PSEPlp_chgobj(&m_lp, ecount, &lp_indices[0], &tour_obj[0]);
  if(rval) return rval;

  rval = PSEPlp_primal_opt(&m_lp, &infeas);
  if(rval) return rval;

  rval = PSEPlp_objval(&m_lp, &objval);
  if(objval != -ncount){
    cerr << "Basis rebuild gave wrong solution\n";
    return rval;
  }

  rval = set_edges();
  if(rval) return rval;

  rval = PSEPlp_getbase(&m_lp, &old_colstat[0], NULL);
  if(rval) return rval;

  prune = PSEP_zeit();
  rval = prune_cuts(&num_removed);
  if(rval) return rval;
  prune = PSEP_zeit() - prune;
  cout << "Removed " << num_removed << " cuts from the LP in "
       << prune << "s\n";


  rval = PSEPlp_chgobj(&m_lp, ecount, &lp_indices[0], &old_obj[0]);
  if(rval) return rval;
    
  return rval;
}

int PSEP_LP_Core::basis_init(){
  old_colstat.resize(m_graph.edge_count);
  for(int i = 0; i < old_colstat.size(); i++)
    old_colstat[i] = CPX_AT_LOWER;
  
  for(int i = 0; i < m_graph.edge_count; i++){
    Edge e = m_graph.edges[i];
    int ind0, ind1;
    if(perm[e.end[0]] < perm[e.end[1]]){
      ind0 = perm[e.end[0]]; ind1 = perm[e.end[1]];
    } else {
      ind1 = perm[e.end[0]]; ind0 = perm[e.end[1]];
    }
      
    if(ind1 - ind0 == 1 || (ind0 == 0 && ind1 == m_graph.node_count - 1)){
      old_colstat[i] = CPX_BASIC;
      if(m_graph.node_count %2 == 1)
	continue;	
    }

    if(m_graph.node_count % 2 == 0){
      if(ind0 == m_graph.node_count - 2 && ind1 == m_graph.node_count - 1){
	old_colstat[i] = CPX_AT_UPPER;
	continue;
      }

      if(ind0 == 0 && ind1 == m_graph.node_count - 2){
	old_colstat[i] = CPX_BASIC;
      }
    } 
  }
  
  int rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);
  if(rval) goto CLEANUP;

  rval = pivot();
  if(rval) goto CLEANUP;
    
  rval = set_edges();
  if(rval) goto CLEANUP;  
  
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
  int num_removed;
  
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

  if(prune_cuts(&num_removed)){
    cerr << "Error entry point: LP_Core::update_best_tour\n";
    return 1;
  } else
    cout << num_removed << " slack cuts pruned from LP after augmenting\n";

  return 0;
};

int PSEP_LP_Core::pivot_until_change(int *pivot_status_p){
int rval = 0;
  int itcount = 0, icount = 0;
  int rowcount = PSEPlp_numrows(&m_lp), ncount = m_graph.node_count;  
  bool integral = false, conn = false, dual_feas = false;
  *pivot_status_p = -1;

  double round_start = PSEP_zeit();

  bool did_jumpstart = false;

  old_rowstat.resize(rowcount);

  while(true){
    if(++itcount == 3 * ncount)
      if(prefs.switching_choice == LP::PRICING::SWITCHING::DYNAMIC)
	change_pricing();

    if(itcount > 3 * ncount){
      if(!did_jumpstart &&
	 prefs.switching_choice == LP::PRICING::SWITCHING::OFF)
	if(prefs.jumpstart){
	  cout << "Temporarily engaging steepest edge...";
	  enable_jumpstart();
	  did_jumpstart = true;
	}
    }

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

  if(did_jumpstart){
    cout << "Now reverting to original pricing\n";
    disable_jumpstart();
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

    // cout << "Did " << itcount << " pivots in "
    // 	 << setprecision(2) << round_start << "s, "
    // 	 << "obj val: " << setprecision(6) << get_obj_val() << "\n";

 CLEANUP:
    if(rval)
      cerr << "Error entry point: pivot_until_change" << endl;
    return rval;
}

void PSEP_LP_Core::change_pricing(){
  int newprice;

  cout << "//////SWITCHED PRICING TO ";
  switch(prefs.pricing_choice){
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

  if(CPXsetintparam(m_lp.cplex_env, CPXPARAM_Simplex_PGradient, newprice))
    cerr << "ERROR: PRICING SWITCH DID NOT TAKE PLACE\n";
  prefs.switching_choice = LP::PRICING::SWITCHING::OFF;
}

void PSEP_LP_Core::enable_jumpstart(){
  if(CPXsetintparam(m_lp.cplex_env, CPXPARAM_Simplex_PGradient,
		    CPX_PPRIIND_STEEP))
    cerr << "ERROR: JUMPSTART PRICING DID NOT TAKE PLACE\n";
}

void PSEP_LP_Core::disable_jumpstart(){
  int oldprice;
  switch(prefs.switching_choice){
  case LP::PRICING::SWITCHING::OFF:
    oldprice = CPX_PPRIIND_PARTIAL;
    break;
  default:
    switch(prefs.pricing_choice){
    case LP::PRICING::DEVEX:
      oldprice = CPX_PPRIIND_DEVEX;
      break;
    case LP::PRICING::STEEPEST:
      oldprice = CPX_PPRIIND_STEEPQSTART;
      break;
    }
  }

  if(CPXsetintparam(m_lp.cplex_env, CPXPARAM_Simplex_PGradient,
		    oldprice))
    cerr << "ERROR: JUMPSTART UNDO DID NOT TAKE PLACE\n";
}

int PSEP_LP_Core::prune_cuts(int *num_removed){
  int rval = 0;
  *num_removed = 0;
  int ncount = best_tour_nodes.size();
  int rowcount = PSEPlp_numrows(&m_lp);
  vector<int> delset(rowcount, 0);
  vector<double> slacks(rowcount - ncount, 0);

  rval = PSEPlp_getslack(&m_lp, &slacks[0], ncount, rowcount - 1);
  if(rval)
    return 1;

  for(int i = 0; i < slacks.size(); i++){
    if(fabs(slacks[i]) >= LP::EPSILON){
      delset[ncount + i] = 1;
      (*num_removed)++;
    }
  }

  rval = PSEPlp_delsetrows(&m_lp, &delset[0]);
  if(rval)
    return 1;

  return 0;
}

