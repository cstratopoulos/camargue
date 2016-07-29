#include "purecut.h"

#include<iomanip>

using namespace std;

int PSEP_PureCut::solve(const bool heur){
  int rval = 0;

  int stat;
  int num_seg, num_2match, num_dp, num_bad, total_cuts = 0;
  int segval, matchval, dpval;
  double segtime, matchtime, dptime, pivtime;
  double piv_val;
  int rounds = 0, augrounds = 0;
  bool in_subtour;

  int max_per_round = LPcore.prefs.max_cuts_round;

  double total_segtime = 0, total_2mtime = 0, total_dptime = 0;
  int total_segcalls = 0, total_2mcalls = 0;
  double total_pivtime = 0;
  double routine_start, fixing_start;
  
  bool called_heur = false;
  bool fixing = LPcore.prefs.redcost_fixing;

  if(fixing){
    fixing_start = PSEP_zeit();
    rval = LPfix.redcost_fixing();
    if(rval) goto CLEANUP;
    fixing_start = PSEP_zeit() - fixing_start; //previously called rebuild
  }
  
  rval = LPcore.basis_init();
  if(rval) goto CLEANUP;
  
  cout << "Pivoting until optimality or no more cuts" << endl;
  routine_start = PSEP_zeit();
  while(/*true*/++rounds < 50){
    //    rounds++;
    augrounds++;
    in_subtour = false;

    pivtime = PSEP_zeit();
    rval = LPcore.pivot_until_change(&stat);
    if(rval) goto CLEANUP;
    total_pivtime += PSEP_zeit() - pivtime;

    if(rounds % 10 == 0)
      piv_val = LPcore.get_obj_val();
    //print.pivot(stat);

    if(stat == PIVOT::FATHOMED_TOUR){
      cout << "\n\n    ROUND " << rounds << " -- ";
      print.pivot(stat);
      if(heur){
	if(Aug.active()){
	  rval = Aug.clear_clamps();
	  if(rval)
	    goto CLEANUP;
	  called_heur = true;
	  if(LPcore.rebuild_basis())
	    goto CLEANUP;
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
	if(Aug.active())
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

    if(rounds % 10 == 0){
      cout << "\n PIVOTING ROUND: " << rounds << " [ "
	   << (LPcore.numrows() - LPcore.best_tour_nodes.size())
	   << " cuts in the LP ]\n";
      cout << "   ";
      print.pivot(stat);
      cout << "   Pivot objval: " << piv_val << "\n";
      cout << "   Avg piv time: " << setprecision(2)
	   << ((double) (total_pivtime / rounds)) <<"s\n"
	   << setprecision(6);
	   
    }

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
	  total_dptime += dptime;
	  cout << "  Round: " << rounds << ", ";
	  cout << "added " << num_dp << " simple DP cuts "
	       << "in " << setprecision(2) << dptime << "s" << setprecision(6);
	  if(num_bad != 0)
	    cout << "!! " << num_bad << " bad cuts found)";
	  cout << "\n";
	}
      }
    }

    if(augrounds >= 25 && (augrounds % 15) == 10
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

  if(stat != 3 || rounds == 50)
    cout << "Terminated due to lack of cutting planes after "
	 << rounds << " rounds of separation" << endl;
  cout << "\n";
  cout << "  " << total_cuts << " cutting planes added over "
       << rounds << " rounds of separation" << endl;
  // cout << "Average time per non-degenerate pivot: "
  //      << ((double) (total_pivtime / rounds)) << "\n";

  cout << "   Total time during lightDP sep: "
       << total_dptime << "s\n";
  cout << "   Average time per segment call: "
       << ((double) (total_segtime / total_segcalls)) << "\n";
  cout << "                     2match call: "
       << ((double) (total_2mtime / total_2mcalls)) << "\n";


  if(called_heur){
    cout << "VVV SOME PRICING OF CLAMPED EDGES SHOULD GO HERE? VVV " << endl;
  }
    
  cout << "\n Total time for Purecut::solve: "
       << (PSEP_zeit() - routine_start) << "\n";
  if(fixing)
    cout <<"         LPfix::redcost_fixing: "
	 << fixing_start << "s\n";

  if(stat != PIVOT::FATHOMED_TOUR)
    LPfix.redcost_fixing();


 CLEANUP:
  if(rval)
    cerr << "Error entry point: PureCut::solve()\n";
  return rval;
}
