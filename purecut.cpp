#include "purecut.h"

using namespace std;

int PSEP_PureCut::solve(const bool heur){
  int rval = 0;

  int stat;
  int num_seg, num_2match, num_dp, num_bad, total_cuts = 0;
  int segval, matchval, dpval;
  double segtime, matchtime, dptime, pivtime;
  int rounds = 0, augrounds = 0;
  bool in_subtour;

  int max_per_round = LPcore.prefs.max_cuts_round;

  double total_segtime = 0, total_2mtime = 0;
  int total_segcalls = 0, total_2mcalls = 0;
  double total_pivtime = 0;

  rval = LPcore.basis_init();
  if(rval) goto CLEANUP;

  cout << "Pivoting until optimality or no more cuts" << endl;

  while(true){
    rounds++;
    augrounds++;
    in_subtour = false;

    pivtime = PSEP_zeit();
    rval = LPcore.pivot_until_change(&stat);
    if(rval) goto CLEANUP;
    total_pivtime += PSEP_zeit() - pivtime;

    print.pivot(stat);

    if(stat == PIVOT::FATHOMED_TOUR){
      if(heur){
	if(Aug.active()){
	  stat = PIVOT::TOUR;

	  cout << "!!!! Treating fathom like augment bc clamps !!!!!!\n";
	  if(LPcore.update_best_tour())
	    goto CLEANUP;
	  
	  rval = Aug.clear_clamps();
	  if(rval)
	    goto CLEANUP;
	  cout << "Now Aug.active() is: " << Aug.active() << "\n";
	  augrounds = 0;
	  continue;
	}
      }
      
      break;
    }

    if(stat == PIVOT::TOUR){
      cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "!!!AUGMENTED TOUR!!!!" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "~Call to delete slack cuts should go here~" << endl;
      if(LPcore.update_best_tour())
	goto CLEANUP;
      if(heur)
	if(Aug.clear_clamps())
	  goto CLEANUP;
      augrounds = 0;
      continue;
    }

    rval = LPcore.pivot_back();
    if(rval) goto CLEANUP;

    segtime = PSEP_zeit();
    segval = cutcall.segment(max_per_round, &num_seg);
    if(segval == 1)
      break;
    segtime = PSEP_zeit() - segtime;
    total_segtime += segtime; total_segcalls++;

    matchtime = PSEP_zeit();
    matchval = cutcall.blossom(max_per_round - num_seg, &num_2match);
    if(matchval == 1)
      break;
    matchtime = PSEP_zeit() - matchtime;
    total_2mtime += matchtime; total_2mcalls++;

    cout << "Added " << num_seg << " segments"
	 << " (in " << segtime << "s)"
	 <<", " << num_2match
	 << " blossoms"
	 << " (in " << matchtime << "s)" << endl;

    dpval = 2; num_dp = 0;
    if(augrounds >= LPcore.prefs.dp_threshold){
      if(num_seg == 0 && stat != PIVOT::SUBTOUR){
	rval = cutcall.in_subtour_poly(&in_subtour);
	if(rval) goto CLEANUP;

	if(in_subtour){
	  dptime = PSEP_zeit();
	  dpval = cutcall.simpleDP(max_per_round - num_2match,
				   &num_dp, &num_bad);
	  if(dpval == 1)
	    break;
	  dptime = PSEP_zeit() - dptime;
	  cout << "*** Added " << num_dp << " simple DP cuts "
	       << "(in " << dptime << "s) (" << num_bad << " bad cuts found)"
	       << " ***\n";
	}
      }
    }

    if(augrounds >= 25 && (augrounds % 25) == 0
       && stat != PIVOT::SUBTOUR){
      if(heur){
	rval = Aug.add_clamp();
	if(rval) goto CLEANUP;
      }
    }

    total_cuts += num_seg + num_2match + num_dp;

    if(num_seg == 0 && num_2match == 0 && num_dp == 0)
      break;
  }
   
  if(stat != 3)
    cout << "Terminated due to lack of cutting planes after "
	 << rounds << " rounds of separation" << endl;
  cout << total_cuts << " cutting planes added over "
       << rounds << " rounds of separation" << endl;
  cout << "Average time per non-degenerate pivot: "
       << ((double) (total_pivtime / rounds)) << "\n";

  cout << "Average time per segcall: "
       << ((double) (total_segtime / total_segcalls)) << "\n";
  cout << "Average time per blossom call: "
       << ((double) (total_2mtime / total_2mcalls)) << "\n";

  rval = LPcore.primal_opt();
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval)
    cerr << "Error entry point: PureCut::solve()\n";
  return rval;
}
