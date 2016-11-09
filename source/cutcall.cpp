#include "cutcall.hpp"
#include "DPgraph.hpp"

#include <chrono>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;


//#define PSEP_TEST_TOOTH

namespace PSEP {

int CutControl::primal_sep(const int augrounds, const LP::PivType stat)
{
  int rval = 0;
  
  int segval = 2, matchval = 2, dpval = 2;
  bool pool_blossoms;

  std::chrono::time_point<std::chrono::system_clock> match_start, match_end;
  std::chrono::duration<double> match_elapsed;

  segtime.resume();
  segval = segments.cutcall();
  if(segval == 1){
    rval = 1;
    goto CLEANUP;
  }  
  segtime.stop();
  
  total_segcalls++;

  matchtime.resume();
  rval = q_has_viol(pool_blossoms, blossom_q);
  if(rval) goto CLEANUP;

  if(pool_blossoms){
    blossom_q.q_fresh = false;
  }
  else {
    matchval = fastblossoms.cutcall();
    if(matchval == 1){ rval = 1; goto CLEANUP; }

    /**@todo old version was matchval == 2. maybe should be a switch? */
    if(matchval == 2 && segval == 2){
      matchval = blossoms.cutcall();
      if(matchval == 1){ rval = 1; goto CLEANUP; }
    }
    
    blossom_q.q_fresh = true;
  }

  matchtime.stop();
  
  total_2mcalls++;

#ifdef PSEP_TEST_TOOTH
  if(segval == 2 && stat != LP::PivType::Subtour){
    bool in_sub = false;
    rval = in_subtour_poly(in_sub);
    if(rval) goto CLEANUP;

    if(in_sub){
      cout << "Augrounds " << augrounds << ", solution is in subtour polytope, "
	   << "testing candidate teeth!" << endl;

      rval = candidates.get_light_teeth();
      if(rval) goto CLEANUP;

      cout << "Got collection of light candidate teeth" << endl;

      DPCutGraph cutgraph(candidates.light_teeth, G_s);

      rval = 1; goto CLEANUP;
    }
  }
#endif

  if(segval == 2 && matchval == 2 && dpval == 2)
    rval = 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in CutControl::primal_sep: ";
  if(segval == 1)
    cerr << "Cuts<seg>::cutcall\n";
  if(matchval == 1)
    cerr << "Cuts<blossom>::cutcall\n";
  if(dpval == 1)
    cerr << "Cuts<domino>::cutcall\n";
  return rval;
}

int CutControl::add_primal_cuts()
{
  int rval = 0;
  int seg_added = 0, blossom_added = 0;
  vector<int> rmatind;
  vector<double> rmatval;
  char sense;
  double rhs;
  int rmatbeg = 0;

  while(!segment_q.empty() && seg_added < prefs.max_per_round){
    rval = translator.get_sparse_row(segment_q.peek_front(),
				     rmatind, rmatval, sense, rhs);
    if(rval) goto CLEANUP;

    rval = PSEPlp_addrows(&m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
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

      rval = PSEPlp_addrows(&m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
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
	rval = PSEPlp_addrows(&m_lp, 1, rmatind.size(), &rhs, &sense, &rmatbeg,
			      &rmatind[0], &rmatval[0]);
	if(rval) goto CLEANUP;

	
	blossom_added++;
      }
      blossom_q.pop_front();
    } 
  }


 CLEANUP:
  if(rval)
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
  dptime.report(false);
  gmitime.report(false);
  cout << "\n";
}

int CutControl::q_has_viol(bool &result, CutQueue<HyperGraph> &pool_q)
{
  int rval = 0;
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
  int ecount = support_ecap.size(), ncount = G_s.node_count;
  double cutval = 2;
  result = false;

  if(CCcut_mincut(ncount, ecount, &support_elist[0], &support_ecap[0],
		  &cutval, NULL, NULL)){
    std::cerr << "Problem in CutControl::in_subtour_poly with CCcut_mincut\n";
    return 1;
  }
    
  if(cutval >= 2)
    result = true;
  
  return 0;
}

}
