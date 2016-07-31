#include "purecut.h"

#include<iomanip>

using namespace std;

int PSEP_PureCut::solve(){
  int rval = 0;

  int stat;
  int num_seg = 0, num_2match = 0, num_dp = 0, total_cuts = 0;
  int segval = 2, matchval = 2, dpval = 2;
  double segtime, matchtime, dptime, pivtime;
  double piv_val;
  int rounds = 0, augrounds = 0;

  double total_segtime = 0, total_2mtime = 0, total_dptime = 0;
  int total_segcalls = 0, total_2mcalls = 0;
  double total_pivtime = 0, max_pivtime = 0;
  double routine_start, fixing_start;

  int roundlimit = 500;

  bool fixing = LPcore.prefs.redcost_fixing;

  if(fixing){
    fixing_start = PSEP_zeit();
    rval = LPfix.redcost_fixing();
    if(rval) goto CLEANUP;
    fixing_start = PSEP_zeit() - fixing_start; 
  }
  
  rval = LPcore.basis_init();
  if(rval) goto CLEANUP;
  
  cout << "Pivoting until optimality or no more cuts" << endl;
  routine_start = PSEP_zeit();
  while(++rounds < roundlimit){
    augrounds++;

    if(rounds % 50 == 0){
      cout << "Calling edge elimination again...\n\n  ";
      rval = LPfix.redcost_fixing();
      if(rval) goto CLEANUP;

      rval = LPcore.rebuild_basis();
      if(rval) goto CLEANUP;            
    }

    pivtime = PSEP_zeit();
    rval = LPcore.pivot_until_change(&stat);
    if(rval) goto CLEANUP;
    pivtime = PSEP_zeit() - pivtime;
    total_pivtime += pivtime;
    if(pivtime > max_pivtime) max_pivtime = pivtime;

    if(rounds % 10 == 0)
      piv_val = LPcore.get_obj_val();
    //print.pivot(stat);

    if(stat == PIVOT::FATHOMED_TOUR){
      cout << "\n\n    ROUND " << rounds << " -- ";
      print.pivot(stat);
      cout << "                Pivot objval: "
	   << LPcore.get_obj_val() << "\n";      
      break;
    }

    if(stat == PIVOT::TOUR){
      cout << "\n\n    !!!AUGMENTED TOUR!!!!" << endl;
      print.pivot(stat);
      cout << "                Pivot objval: "
	   << LPcore.get_obj_val() << "\n"; 
      if(LPcore.update_best_tour())
	goto CLEANUP;
      augrounds = 0;
      continue;
    }

    rval = LPcore.pivot_back();
    if(rval) goto CLEANUP;

    num_seg = 0;
    segtime = PSEP_zeit();
    segval = CutControl.segments.cutcall();
    if(segval == 1){
      rval = 1;
      goto CLEANUP;
    }
    if(segval == 0)
      num_seg = 1;
    segtime = PSEP_zeit() - segtime;
    total_segtime += segtime; total_segcalls++;

    num_2match = 0;
    matchtime = PSEP_zeit();
    matchval = CutControl.blossoms.cutcall();
    if(matchval == 1){
      rval = 1;
      goto CLEANUP;
    }
    if(matchval == 0) num_2match = 1;
    matchtime = PSEP_zeit() - matchtime;
    total_2mtime += matchtime; total_2mcalls++;

    dpval = 2; num_dp = 0;
    if(augrounds >= LPcore.prefs.dp_threshold){
      if(num_seg == 0 && stat != PIVOT::SUBTOUR){
	dptime = PSEP_zeit();
	dpval = CutControl.dominos.cutcall();
	if(dpval == 1){
	  rval = 1; goto CLEANUP;
	}
	dptime = PSEP_zeit() - dptime;
	total_dptime += dptime;
      }	 
    }


    if(rounds % 10 == 0){
      cout << "\n PIVOTING ROUND: " << rounds << " [ "
	   << (LPcore.numrows() - LPcore.best_tour_nodes.size())
	   << " cuts in the LP ]\n";
      cout << "   ";
      print.pivot(stat);
      cout << "   Pivot objval: " << piv_val << "\n";
      cout << "   Avg piv time: " << setprecision(2)
	   << ((double) (total_pivtime / rounds)) <<"s, (longest "
	   << max_pivtime << "s)\n" << setprecision(6);
      max_pivtime = 0;
	   
    }

    total_cuts += num_seg + num_2match + num_dp;

    if(num_seg == 0 && num_2match == 0 && num_dp == 0)
      break;
  }

  if(stat != 3 || rounds == roundlimit){
    cout << "\n Terminated due to: ";
    if(num_seg == 0 && num_2match == 0 && num_dp == 0){
      cout << "lack of cutting planes.\n      ";
      print.pivot(stat);
      cout << "\n";
    }
    else if(rounds == roundlimit)
      cout << "artificial round limit\n";
  }
  
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

    
  cout << "\n Total time for Purecut::solve: "
       << (PSEP_zeit() - routine_start) << "s\n";
  if(fixing)
    cout <<"         LPfix::redcost_fixing: "
	 << fixing_start << "s\n";

  if(stat != PIVOT::FATHOMED_TOUR)
    LPfix.redcost_fixing();


 CLEANUP:
  if(segval == 1)
    cerr << "Problem in Cuts<seg>::cutcall()\n";
  if(matchval == 1)
    cerr << "Problem in Cuts<blossom>::cutcall()\n";
  if(rval)
    cerr << "Error entry point: PureCut::solve()\n";
  return rval;
}
