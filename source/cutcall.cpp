#include "cutcall.hpp"
#include "DPgraph.hpp"
#include "util.hpp"

#include <iomanip>
#include <algorithm>

extern "C" {
 #include <concorde/INCLUDE/combs.h>
}

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

namespace CMR {

int CutControl::primal_sep(const int augrounds, const LP::PivType stat)
{
  int rval = 0;
  
  int segval = 2, matchval = 2, dpval = 2;
  bool pool_blossoms;

  segtime.resume();
  segval = segments.cutcall();
  if(segval == 1){ rval = 1; goto CLEANUP; }  
  segtime.stop();
  
  total_segcalls++;

  matchtime.resume();
  rval = q_has_viol(pool_blossoms, blossom_q);
  if(rval) goto CLEANUP;

  if(pool_blossoms){
    blossom_q.q_fresh = false;
  }
  else {
    /**@todo old version was matchval == 2. maybe should be a switch? */
    if(matchval == 2 && segval == 2){
      matchval = blossoms.cutcall();
      if(matchval == 1){ rval = 1; goto CLEANUP; }
    }
    
    blossom_q.q_fresh = true;
  }

  matchtime.stop();
  
  total_2mcalls++;

  /*
  if(segval == 2 && matchval == 2 && stat != LP::PivType::Subtour){
    bool in_sub = false;
    
    rval = in_subtour_poly(in_sub);
    if(rval) goto CLEANUP;

    if(in_sub){
      dptime.resume();
      try {
	dominos = CMR::make_unique<Cut<dominoparity>>(graph_data,best_data,
						       supp_data, domino_q);
      } catch(...){ CMR_SET_GOTO(rval, "Couldn't allocate dominos. "); }

      dpval = dominos->cutcall();
      if(dpval == 1){ rval = 1; goto CLEANUP; }
      dptime.stop();

      if(dpval == 0){
	while(!domino_q.empty()){
	  vector<int> rmatind;
	  vector<double> rmatval;
	  char sense;
	  double rhs, tour_activity;

	  rval = translator.get_sparse_row(domino_q.peek_front(),
					   best_data.best_tour_nodes, rmatind,
					   rmatval, sense, rhs);
	  if(rval) goto CLEANUP;

	  translator.get_activity(tour_activity, best_data.best_tour_edges,
				  rmatind, rmatval);
	  if(tour_activity == rhs)
	    break;
	  else
	    domino_q.pop_front();
	}

	if(domino_q.empty()) dpval = 2;
      }
    }
  }
  */

  if(segval == 2 && matchval == 2 && dpval == 2)
    rval = 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in CutControl::primal_sep: ";
  if(segval == 1)
    cerr << "Cuts<seg>::cutcall";
  if(matchval == 1)
    cerr << "Cuts<blossom>::cutcall";
  if(dpval == 1)
    cerr << "Cuts<domino>::cutcall";
  if(rval == 1)
    cerr << "\n";  
  return rval;
}

int CutControl::add_primal_cuts()
{
  int rval = 0;
  int seg_added = 0, blossom_added = 0, dp_added = 0;
  vector<int> rmatind;
  vector<double> rmatval;
  char sense;
  double rhs;
  int rmatbeg = 0;

  LP::Prefs &prefs = LP_data.prefs;
  CMRlp *m_lp = &LP_data.m_lp;
  vector<double> &m_lp_edges = LP_data.m_lp_edges;

  while(!segment_q.empty() && seg_added < prefs.max_per_round){
    rval = translator.get_sparse_row(segment_q.peek_front(),
				     rmatind, rmatval, sense, rhs);
    if(rval) goto CLEANUP;

    rval = CMRlp_addrows(m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
			  &rmatind[0], &rmatval[0]);
    if(rval) goto CLEANUP;

    seg_added++;
    segment_q.pop_front();
  }

  if(blossom_q.q_fresh){
    while(!blossom_q.empty() && blossom_added < prefs.max_per_round){
      rval = translator.get_sparse_row(blossom_q.peek_front(),
				       rmatind, rmatval, sense, rhs);
      if(rval) goto CLEANUP;

      rval = CMRlp_addrows(m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
			    &rmatind[0], &rmatval[0]);
      if(rval) goto CLEANUP;

      blossom_added++;
      blossom_q.pop_front();
    } 
  } else {
    while(!blossom_q.empty() && blossom_added < prefs.max_per_round){
      bool is_violated = false;
      rval = translator.get_sparse_row_if(is_violated, blossom_q.peek_front(),
					  m_lp_edges, rmatind, rmatval,
					  sense, rhs);
      if(rval) goto CLEANUP;

      if(is_violated){
	rval = CMRlp_addrows(m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
			      &rmatind[0], &rmatval[0]);
	if(rval) goto CLEANUP;

	
	blossom_added++;
      }
      blossom_q.pop_front();
    } 
  }

  while(!domino_q.empty()){
    const dominoparity &dp_cut = domino_q.peek_front();
    double tour_activity;

    rval = translator.get_sparse_row(dp_cut, best_data.best_tour_nodes,
				     rmatind, rmatval, sense, rhs);
    if(rval) goto CLEANUP;

    translator.get_activity(tour_activity, best_data.best_tour_edges,
			    rmatind, rmatval);

    if(tour_activity == rhs){
      rval = CMRlp_addrows(m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
			    &rmatind[0], &rmatval[0]);
      if(rval) goto CLEANUP;
      
      ++dp_added;
    }

    domino_q.pop_front();
  }

  if(dp_added)
    cout << "\tAdded " << dp_added << " simple DP inequalities\n";

 CLEANUP:
  if(rval == 1)
    cerr << "CutControl::add_primal_cuts failed\n";
  return rval;
}

int CutControl::safe_gomory_sep(){
  gmitime.resume();
  int rval = safe_gomory.cutcall();
  gmitime.stop();

  total_gencalls++;
  if(rval == 1)
    cerr << "Problem in CutControl::safe_gomory_sep\n";
  return rval;
}

void CutControl::profile()
{
  segtime.report(false);
  matchtime.report(true);
  dptime.report(true);
  gmitime.report(false);
  cout << "\n";
}

int CutControl::q_has_viol(bool &result, CutQueue<HyperGraph> &pool_q)
{
  int rval = 0;
  vector<double> &m_lp_edges = LP_data.m_lp_edges;
  result = false;

  while(!pool_q.empty()){
    rval = translator.is_cut_violated(result, pool_q.peek_front(),
				      m_lp_edges);
    if(rval) goto CLEANUP;

    if(result)
      break;

    pool_q.pop_front();
  }

 CLEANUP:
  if(rval)
    cerr << "CutControl::q_contains_viol failed\n";
  return rval;
}

int CutControl::in_subtour_poly(bool &result)
{
  vector<int> &support_elist = supp_data.support_elist;
  vector<double> &support_ecap = supp_data.support_ecap;

  int ecount = support_ecap.size(), ncount = supp_data.G_s.node_count;
  double cutval = 2;
  double rhs = 2.0 - Epsilon::Cut;

  result = false;

  if(CCcut_mincut(ncount, ecount, &support_elist[0], &support_ecap[0],
		  &cutval, nullptr, nullptr)){
    cerr << "Problem in CutControl::in_subtour_poly with CCcut_mincut.\n";
    return 1;
  }

  

  result = (cutval > rhs);

  return 0;
}

}
