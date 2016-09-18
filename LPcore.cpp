#include "LPcore.h"

#include <math.h>
#include <iomanip>

using namespace std;
using namespace PSEP::LP;

bool Core::is_dual_feas(){
  return PSEPlp_dualfeas(&m_lp);
}

inline bool Core::is_integral(){
  for(int i = 0; i < m_lp_edges.size(); i++)
    if((m_lp_edges[i] >= EPSILON) && (m_lp_edges[i] <= (1 - EPSILON)))
      return false;

  return true;
}

int Core::is_best_tour_feas(bool &result){
  int rval = 0;
  vector<double> feas_stat;
  result = true;

  if(best_tour_edges_lp.size() != best_tour_edges.size()){
      best_tour_edges_lp.resize(best_tour_edges.size());
      for(const int i: best_tour_edges)
	best_tour_edges_lp[i] = best_tour_edges[i];
    }

  try { feas_stat.resize(numrows()); } catch(const std::bad_alloc &) {
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for feas stat, ");
  }

  rval = PSEPlp_getrowinfeas(&m_lp, &best_tour_edges_lp[0], &feas_stat[0],
			     0, numrows() - 1);
  if(rval) goto CLEANUP;

  for(int i = 0; i < feas_stat.size(); i++)
    if(feas_stat[i] != 0.0){
      cout << "Row " << i << " infeasible, feas_stat: "
	   << feas_stat[i] << "\n";
      result = false;
    }

 CLEANUP:
  if(rval)
    cerr << "Problem in Core::is_best_tour_feas\n";
  return rval;
}

int Core::factor_basis(){
  int rval = PSEPlp_no_opt(&m_lp);
  if(rval)
    cerr << "Problem in LPCore::factor_basis\n";
  return rval;
}

int Core::single_pivot(){
  int infeasible = 0;
  int rval = PSEPlp_primal_pivot(&m_lp, &infeasible);

  if(rval || infeasible)
    cerr << "Problem in LP_Core::single_pivot(), infeasible "
	 << infeasible << "\n";
  return rval;
}

int Core::nondegenerate_pivot(){
  int infeasible = 0, rval = 0;
  double lowlimit = m_min_tour_value - EPSILON;

  rval = PSEPlp_primal_nd_pivot(&m_lp, &infeasible, lowlimit);
  if(rval || infeasible){
    cerr << "Problem in LPCore::nondegenerate_pivot, infeasible: "
	 << infeasible << "\n";
    rval = 1;
  }
  return rval;
}

int Core::primal_opt(){
  int infeasible = 0;
  double start = zeit();
  int rval = PSEPlp_primal_opt(&m_lp, &infeasible);
  start = zeit() - start;

  if(rval)
    cerr << "Entry point LP_Core::primal_opt(), infeasible "
	 << infeasible << "\n";
  else
    cout << "Primal optimized with obj val " << get_obj_val()
	 << " in " << start << "s.\n";

  return rval;
}


int Core::pivot_back(){
  bool tour_feas = false;
  //int rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);
  int rval = PSEPlp_copystart(&m_lp, &old_colstat[0], &old_rowstat[0],
			      &best_tour_edges_lp[0], NULL, NULL, NULL);
  if(rval) goto CLEANUP;

  rval = factor_basis();
  if(rval) goto CLEANUP;

  rval = is_best_tour_feas(tour_feas);
  if(rval) goto CLEANUP;
  
  if(!tour_feas){
    cout << "Best tour is infeasible after pivot back!\n";
    rval = 1;
  }

 CLEANUP:
  if(rval)
    cerr << "Error entry point: pivot_back()" << endl;
  return rval;
}

double Core::get_obj_val(){
  double objval;
  PSEPlp_objval(&m_lp, &objval);
  return objval;
}

double Core::set_edges(){
  int rval = 0;

  if (m_lp_edges.size() != m_graph.edge_count)
    m_lp_edges.resize(m_graph.edge_count);

  rval = PSEPlp_x (&m_lp, &m_lp_edges[0]);
  if(rval)
    fprintf(stderr, "failed to set_edges(), rval %d\n", rval);

  return rval;
}

int Core::rebuild_basis(bool prune){
  int rval = 0;
  int ecount = m_lp_edges.size();
  int num_removed = 0;
  bool tour_feas;
  
  old_colstat.resize(ecount);
  best_tour_edges_lp.resize(ecount);

  //  double rebuild_time = zeit();
  for(int i = 0; i < m_lp_edges.size(); i++)
    best_tour_edges_lp[i] = best_tour_edges[i];

  rval = PSEPlp_copystart(&m_lp,
			  NULL, NULL,
			  &m_lp_edges[0], NULL,
			  NULL, NULL);
  if(rval) goto CLEANUP;

  rval = factor_basis();
  if(rval) goto CLEANUP;
  
  rval = set_edges();
  if(rval) goto CLEANUP;

  tour_feas = false;
  rval = is_best_tour_feas(tour_feas);
  if(rval) goto CLEANUP;

  if(tour_feas)
    cout << "Best tour is still feasible\n";
  else {
    cout << "Best tour is infeasible!\n";
    rval = 1; goto CLEANUP;
  }


  rval = PSEPlp_getbase(&m_lp, &old_colstat[0], NULL);
  if(rval) goto CLEANUP;

  if(prune){
    rval = LPPrune.prune_cuts(num_removed);
    if(rval) goto CLEANUP;
    cout << "Removed " << num_removed << " cuts from the LP\n";
  }
    
 CLEANUP:
  if(rval)
    cerr << "Error entry point: LPCore::rebuild_basis\n";
  return rval;
}

int Core::rebuild_basis(int &numremoved, IntPair skiprange,
				std::vector<int> &delset){
  int rval = 0, infeas = 0;
  int ecount = m_lp_edges.size();
  int ncount = best_tour_nodes.size();
  int num_removed = 0;
  double objval;
  vector<double> tour_obj(ecount);
  vector<double> old_obj(ecount);
  vector<int> lp_indices(ecount);
  old_colstat.resize(ecount);

  rval = PSEPlp_getobj(&m_lp, &old_obj[0], ecount);
  if(rval) goto CLEANUP;

  for(int i = 0; i < ecount; i++){
    tour_obj[i] = 1 - (2 * best_tour_edges[i]);
    lp_indices[i] = i;
  }

  rval = PSEPlp_chgobj(&m_lp, ecount, &lp_indices[0], &tour_obj[0]);
  if(rval) goto CLEANUP;

  rval = PSEPlp_primal_opt(&m_lp, &infeas);
  if(rval) goto CLEANUP;

  rval = PSEPlp_objval(&m_lp, &objval);
  if(objval != -ncount){
    cerr << "Basis rebuild gave wrong solution, objval: "
	 << objval << "\n";
    rval = 1;
    goto CLEANUP;
  }

  rval = set_edges();
  if(rval) goto CLEANUP;

  rval = PSEPlp_getbase(&m_lp, &old_colstat[0], NULL);
  if(rval) goto CLEANUP;


  rval = LPPrune.prune_with_skip(num_removed, skiprange, delset);
  if(rval) goto CLEANUP;



  rval = PSEPlp_chgobj(&m_lp, ecount, &lp_indices[0], &old_obj[0]);
  if(rval) goto CLEANUP;
    
 CLEANUP:
  if(rval)
    cerr << "Error entry point: LPCore::rebuild_basis\n";
  return rval;
}

int Core::basis_init(){
  int rval = 0;
  int ecount = best_tour_edges.size();

  //for now can assume size of old colstat is non-increasing, hence no
  //need for try - catch
  old_colstat.resize(ecount); 
  for(const int i: old_colstat)
    old_colstat[i] = CPX_AT_LOWER;

  try { best_tour_edges_lp.resize(best_tour_edges.size()); }
  catch(const bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for lp best tour edge, ");
  }
  
  for(int i = 0; i < ecount; i++)
    best_tour_edges_lp[i] = best_tour_edges[i];
  
  //the routine for constructing a starting basis, as per Padberg-Hong
  for(int i = 0; i < ecount; i++){
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

  //rval = PSEPlp_copybase(&m_lp, &old_colstat[0], &old_rowstat[0]);
  rval = PSEPlp_copystart(&m_lp, &old_colstat[0], &old_rowstat[0],
			      &best_tour_edges_lp[0], NULL, NULL, NULL);
  if(rval) goto CLEANUP;

  rval = factor_basis();
  if(rval) goto CLEANUP;


  { bool result;
  rval = is_best_tour_feas(result);
  if(rval) goto CLEANUP;
  if(result)
    cout << "Best tour feasible\n";
  else {
    cout << "BEST TOUR INFEASIBLE AT BASIS INIT!!!\n";
    rval = 1;
  }
  }
    
  CLEANUP:
  if(rval)
    cerr << "Error entry point: basis_init()" << endl;
  return rval;
}

double Core::set_support_graph(){
  int rval = 0;
  
  support_indices.clear();
  support_elist.clear();
  support_ecap.clear();

  for(int i = 0; i < m_graph.edge_count; i++){
    if(m_lp_edges[i] >= EPSILON){
      support_indices.push_back(i);
      support_ecap.push_back(m_lp_edges[i]);
      support_elist.push_back(m_graph.edges[i].end[0]);
      support_elist.push_back(m_graph.edges[i].end[1]);
    }
  }

  rval = GraphUtils::build_s_graph(m_graph.node_count, support_indices.size(),
				m_graph.edges, support_indices, m_lp_edges,
				&G_s);
  
  if(rval)
    cerr << "Problem in LPCore::set_support_graph" << endl;
  return rval;
}

bool Core::test_new_tour(){
  double objval = 0;
  bool result = false;

  for(int i = 0; i < m_graph.edge_count; i++)
    if(m_lp_edges[i] >= EPSILON)
      objval += m_graph.edges[i].len;

  result = objval < m_min_tour_value && is_integral();
  if(result)
    cout << "Supposed improvement?? objval: " << objval << "\n";
  return result;
}

int Core::update_best_tour(){
  double objval = 0;
  
  for(int i = 0; i < m_graph.node_count; i++)
    best_tour_nodes[i] = island[i];

  for(int i = 0; i < m_graph.edge_count; i++)
    if(m_lp_edges[i] < EPSILON)
      best_tour_edges[i] = 0;
    else {
      best_tour_edges[i] = 1;
      objval += m_graph.edges[i].len;
    }

  if(objval > m_min_tour_value){
    cerr << "New best tour is worse! objval " << objval << "\n";
    return 1;
  }

  m_min_tour_value = objval;

  for(int i = 0; i < m_graph.node_count; i++)
    perm[best_tour_nodes[i]] = i;

  return 0;
}

int Core::pivot_until_change(PivType &pivot_status){
  int rval = 0;
  int icount = 0;
  int rowcount = PSEPlp_numrows(&m_lp);
  bool integral = false, conn = false, dual_feas = false;

  double round_start = zeit();


  old_rowstat.resize(rowcount);

  while(true){
    if((dual_feas = is_dual_feas()))
      break;

    rval = PSEPlp_getbase(&m_lp, &old_colstat[0], &old_rowstat[0]);
    if(rval) goto CLEANUP;

    rval = nondegenerate_pivot();
    if(rval) goto CLEANUP;

    rval = set_edges();
    if(rval) goto CLEANUP;

    if(fabs(get_obj_val() - m_min_tour_value) >= EPSILON)
      break;    
  }
  
  round_start = zeit() - round_start;

  rval = set_support_graph();
  if(rval) goto CLEANUP;

  integral = is_integral();
    if(integral){
    conn = GraphUtils::connected(&G_s, &icount, island, 0);
    if(integral && conn){
      if(dual_feas)
	pivot_status = PivType::FathomedTour;
      else
	pivot_status = PivType::Tour;
    } else
      pivot_status = PivType::Subtour;
    } else{
      pivot_status = PivType::Frac;
      frac_colstat.resize(PSEPlp_numcols(&m_lp));
      frac_rowstat.resize(rowcount);
      rval = PSEPlp_getbase(&m_lp, &frac_colstat[0], &frac_rowstat[0]);
    }

 CLEANUP:
    if(rval)
      cerr << "Error entry point: pivot_until_change" << endl;
    return rval;
}

void Core::change_pricing(){
  int newprice;

  switch(prefs.price_method){
  case Pricing::Devex:
    newprice = CPX_PPRIIND_DEVEX;
    break;
  case Pricing::SlackSteepest:
    newprice = CPX_PPRIIND_STEEPQSTART;
    break;
  case Pricing::Steepest:
    newprice = CPX_PPRIIND_STEEP;
    break;
  }

  if(CPXsetintparam(m_lp.cplex_env, CPXPARAM_Simplex_PGradient, newprice))
    cerr << "ERROR: PRICING SWITCH DID NOT TAKE PLACE\n";
}

