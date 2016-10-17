#include "LPcore.hpp"

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
      for(int i = 0; i < best_tour_edges.size(); i++)
	best_tour_edges_lp[i] = best_tour_edges[i];
    }

  try { feas_stat.resize(numrows()); } catch(const std::bad_alloc &) {
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for feas stat, ");
  }

  rval = PSEPlp_getrowinfeas(&m_lp, &best_tour_edges_lp[0], &feas_stat[0],
			     0, numrows() - 1);
  if(rval) goto CLEANUP;

  for(int i = 0; i < feas_stat.size(); i++)
    //    if(feas_stat[i] != 0.0){
    if(fabs(feas_stat[i]) >= EPSILON){
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
  int rval = PSEPlp_copystart(&m_lp, &old_colstat[0], &old_rowstat[0],
			      &best_tour_edges_lp[0], NULL, NULL, NULL);
  if(rval) goto CLEANUP;

  rval = factor_basis();
  if(rval) goto CLEANUP;

  rval = is_best_tour_feas(tour_feas);
  if(rval) goto CLEANUP;
  
  if(!tour_feas){
    cerr << "Best tour is infeasible after pivot back! "
	 << numrows() << " rows in the LP\n";
    cerr << "Pivot back objval: " << get_obj_val() << endl;
    rval = 1;
  }

 CLEANUP:
  if(rval)
    cerr << "Error entry point: pivot_back()" << endl;
  return rval;
}

int Core::dual_pivot()
{
  int infeasible = 0, rval = 0;
  cout << "PRIMAL PIVOTING INSTEAD LOL!!" << endl;
  rval = nondegenerate_pivot();
  //  rval = PSEPlp_dual_pivot(&m_lp, &infeasible);

  if(rval || infeasible){
    cerr << "Problem in LP::Core::dual_pivot(), infeasible "
	 << infeasible << "\n";
    goto CLEANUP;
  }

  try {
  frac_colstat.resize(numcols());
  frac_rowstat.resize(numrows());
  } catch(...) {
    rval = 1; PSEP_GOTO_CLEANUP("Couldn't resize row/colstat, ");
  }

  rval = PSEPlp_getbase(&m_lp, &frac_colstat[0], &frac_rowstat[0]);
  PSEP_CHECK_RVAL(rval, "Couldn't get frac basis, ");

  rval = set_edges();
  if(rval) goto CLEANUP;

  if(is_integral()){
    rval = 1; PSEP_GOTO_CLEANUP("LP solution still integral! ");
  } else {
    cout << "    Dual pivoted to fractional solution, objval : "
	 << get_obj_val() << ", infeasible: " << infeasible << endl;
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in LP::Core::dual_pivot\n";
  return rval;
}

int Core::add_connect_cuts(PivType &piv_stat)
{
  int rval = 0, deltacount = 0;
  int rmatbeg = 0, newrows = 1;
  char sense = 'G';
  double rhs = 2.0;
  vector<int> island_vec;
  vector<double> rmatval;

  connect_cut_range.first = numrows();
  connect_cut_range.second = connect_cut_range.first;


  do {
    island_vec = vector<int>(island.begin(), island.begin() + icount);
    
    GraphUtils::get_delta(island_vec, m_graph.edges, &deltacount, delta,
			  edge_marks);

    try { rmatval.resize(deltacount, 1.0); } catch(...) {
      rval = 1; PSEP_GOTO_CLEANUP("Couldn't resize rmatval, ");
    }

    connect_cut_range.second = numrows();

    rval = PSEPlp_addrows(&m_lp, newrows, deltacount, &rhs, &sense, &rmatbeg,
			  &delta[0], &rmatval[0]);
    PSEP_CHECK_RVAL(rval, "Couldn't add subtour row, ");

    rval = nondegenerate_pivot();
    if(rval) goto CLEANUP;

    rval = set_edges();
    if(rval) goto CLEANUP;

    rval = set_support_graph();
    if(rval) goto CLEANUP;  
  
  } while(!GraphUtils::connected(&G_s, &icount, island, 0));

  cout << "\n    Added "
       << (connect_cut_range.second - connect_cut_range.first + 1)
       << " connect cuts...";

  if(is_integral()){
    piv_stat = PivType::Tour;
    cout << "\n    !!!Pivoted to new tour!!! ";
  } else {
    piv_stat = PivType::Frac;
    cout << "pivoted to connected fractional solution.";
    try {
      frac_colstat.resize(numcols());
      frac_rowstat.resize(numrows());
    } catch(...) {
      rval = 1; PSEP_GOTO_CLEANUP("Couldn't resize frac stats. ");
    }
    
    rval = PSEPlp_getbase(&m_lp, &frac_colstat[0], &frac_rowstat[0]);
    if(rval) goto CLEANUP;
  }



 CLEANUP:
  if(rval){
    cerr << "Problem in LP::Core::add_connect_cut\n";
    connect_cut_range = IntPair(-1, -1);
  }
  return rval;
}

int Core::del_connect_cut()
{
  int rval = 0;
  if(connect_cut_range.first < 0){
    rval = 1;
    PSEP_GOTO_CLEANUP("Tried to delete negative row, ");
  }

  rval = PSEPlp_delrows(&m_lp, connect_cut_range.first,
			connect_cut_range.second);
  PSEP_CHECK_RVAL(rval, "Couldn't delete subtour row, ");
  cout << "    ...Deleted non-tight connect cut(s)." << endl;

 CLEANUP:
  if(rval)
    cerr << "LP::Core::del_connect_cut failed\n";
  connect_cut_range = IntPair(-1, -1);
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
  double objval;
  
  old_colstat.resize(ecount);
  best_tour_edges_lp.resize(ecount);

  //  double rebuild_time = zeit();
  for(int i = 0; i < m_lp_edges.size(); i++)
    best_tour_edges_lp[i] = best_tour_edges[i];

  rval = PSEPlp_copystart(&m_lp,
			  NULL, NULL,
			  &best_tour_edges_lp[0], NULL,
			  NULL, NULL);
  if(rval) goto CLEANUP;

  rval = factor_basis();
  if(rval) goto CLEANUP;
  
  rval = set_edges();
  if(rval) goto CLEANUP;

  objval = get_obj_val();
  if(fabs(objval - m_min_tour_value) >= EPSILON){
    rval = 1; PSEP_GOTO_CLEANUP("Rebuilt basis and got objval " << objval
				<< ". ");
  }

  tour_feas = false;
  rval = is_best_tour_feas(tour_feas);
  if(rval) goto CLEANUP;

  if(!tour_feas){
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
  int ncount = m_graph.node_count;

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

  for(int i = 0; i < ncount; i++){
    int end0 = fmin(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
    int end1 = fmax(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
    IntPairMap::const_iterator edge_it =
      m_graph.edge_lookup.find(IntPair(end0, end1));
    if(edge_it == m_graph.edge_lookup.end()){
      rval = 1;
      PSEP_GOTO_CLEANUP("Couldn't find " << end0 << ", " << end1 << ", ");
    }

    int edge_index = edge_it->second;
    old_colstat[edge_index] = CPX_BASIC;
  }

  if((ncount % 2) == 0){
    int end0 = fmin(best_tour_nodes[ncount - 2], best_tour_nodes[ncount - 1]);
    int end1 = fmax(best_tour_nodes[ncount - 2], best_tour_nodes[ncount - 1]);

    IntPairMap::const_iterator edge_it =
      m_graph.edge_lookup.find(IntPair(end0, end1));
    if(edge_it == m_graph.edge_lookup.end()){
      rval = 1;
      PSEP_GOTO_CLEANUP("Couldn't find " << end0 << ", " << end1 << ", ");
    }

    //the discarded column, i_{n-1}, i_{n} is at upper
    int edge_index = edge_it->second;
    old_colstat[edge_index] = CPX_AT_UPPER;

    end0 = fmin(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
    end1 = fmax(best_tour_nodes[0], best_tour_nodes[ncount - 2]);

    edge_it = m_graph.edge_lookup.find(IntPair(end0, end1));
    if(edge_it == m_graph.edge_lookup.end()){
      rval = 1;
      PSEP_GOTO_CLEANUP("Couldn't find " << end0 << ", " << end1 << ", ");
    }

    //the edge i_{1}, i_{n-1} is added
    edge_index = edge_it->second;
    old_colstat[edge_index] = CPX_BASIC;
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
  if(!result) {
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

  return result;
}

int Core::update_best_tour(){
  double objval = 0;
  int rval = 0;
  bool newbest_feas = false;
  
  for(int i = 0; i < m_graph.node_count; i++)
    best_tour_nodes[i] = island[i];

  for(int i = 0; i < m_graph.edge_count; i++)
    if(m_lp_edges[i] < EPSILON)
      best_tour_edges[i] = 0;
    else {
      best_tour_edges[i] = 1;
      objval += m_graph.edges[i].len;
    }

  rval = (objval > m_min_tour_value);
  PSEP_CHECK_RVAL(rval, "New best tour is worse! Objval " << objval << "\n");

  m_min_tour_value = objval;

  for(int i = 0; i < m_graph.node_count; i++)
    perm[best_tour_nodes[i]] = i;

  old_rowstat.resize(numrows());

  rval = PSEPlp_getbase(&m_lp, &old_colstat[0], &old_rowstat[0]);
  if(rval) goto CLEANUP;

  rval = is_best_tour_feas(newbest_feas);
  if(rval) goto CLEANUP;

  rval = !newbest_feas;
  if(rval) goto CLEANUP;

  cout << "    Verified new tour is a feasible augmentation\n";

 CLEANUP:
  if(rval)
    cerr << "LP::Core::update_best_tour failed\n";
  return rval;
}

int Core::pivot_until_change(PivType &pivot_status){
  int rval = 0;
  icount = 0;
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
      else {
	pivot_status = PivType::Tour;
      }
    } else {
      pivot_status = PivType::Subtour;
    }
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

