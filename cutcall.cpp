#include "cutcall.hpp"

using namespace std;
using namespace PSEP;

int CutControl::primal_sep(const int augrounds, const LP::PivType stat)
{
  int rval = 0;
  
  int segval = 2, matchval = 2, dpval = 2;
  double segtime, matchtime;//, dptime;
  bool pool_blossoms;

  segtime = zeit();
  segval = segments.cutcall();
  if(segval == 1){
    rval = 1;
    goto CLEANUP;
  }  
  segtime = zeit() - segtime;
  total_segtime += segtime;
  total_segcalls++;


  matchtime = zeit();

  rval = q_has_viol(pool_blossoms, blossom_q);
  if(rval) goto CLEANUP;

  if(pool_blossoms){
    blossom_q.q_fresh = false;
  }
  else {
    matchval = fastblossoms.cutcall();
    if(matchval == 1){ rval = 1; goto CLEANUP; }

    if(matchval == 2){
      matchval = blossoms.cutcall();
      if(matchval == 1){ rval = 1; goto CLEANUP; }
    }
    
    blossom_q.q_fresh = true;
  }

  matchtime = zeit() - matchtime;
  total_2mtime += matchtime;
  total_2mcalls++;

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
      //      candidates.print_collection();

      rval = 1; goto CLEANUP;
    }
  }

  // if(prefs.dp_threshold >= 0 && augrounds > 0 &&
  //    (augrounds % prefs.dp_threshold == 0)){
  //   if(segval == 2 && matchval == 2 && stat != LP::PivType::Subtour){
  //     dptime = zeit();
  //     dpval = dominos.cutcall();
  //     if(dpval == 1){
  // 	rval = 1;
  // 	goto CLEANUP;
  //     }
  //     dptime = zeit() - dptime;
  //     total_dptime += dptime;	
  //   }
  // }

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
  double gentime = zeit();
  int rval = safe_gomory.cutcall();
  gentime = zeit() - gentime;
  total_gentime += gentime;
  total_gencalls++;
  if(rval == 1)
    cerr << "Problem in CutControl::safe_gomory_sep\n";
  return rval;
}

void CutControl::profile(const double total_time)
{
  std::cout << "   Total time during lightDP sep: " << std::setprecision(4)
	    << total_dptime << ", ratio: "
	    << (total_dptime / total_time) << "\n"
	    << "                     segment sep: "
	    << total_segtime << ", ratio: "
	    << (total_segtime / total_time) << "\n"
	    << "                     blossom sep: "
	    << total_2mtime << ", ratio: "
	    << (total_2mtime / total_time) << "\n"		 
	    << "                    safe GMI sep: "
	    << total_gentime << ", ratio: "
	    << (total_gentime / total_time) << "\n"
	    << std::setprecision(6);
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
  int end0 = 0;
  double cutval = 2;
  result = false;

  for(int end1 = 1; end1 < ncount; end1++){
    if(CCcut_mincut_st(ncount, ecount, &support_elist[0], &support_ecap[0],
		       end0, end1, &cutval, (int **) NULL, (int *) NULL)){
      cerr << "CutControl::in_subtour_poly call to Concorde st-cut failed\n";
      return 1;
    }

    if(cutval < 2)
      return 0;
  }

  result = true;
  return 0;
}
