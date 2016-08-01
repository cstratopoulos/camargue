#include "purecut.h"

#include<iomanip>

using namespace std;
using namespace PSEP;

int PureCut::solve(){
  int rval = 0, cut_rval;

  int stat;
  double piv_val;
  int rounds = 0, augrounds = 0;

  double pivtime, total_pivtime = 0, max_pivtime = 0;
  int num_removed = 0;
  double routine_start, fixing_start;

  int roundlimit = 100;

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
      if(LPPrune.prune_cuts(num_removed))
	goto CLEANUP;
      cout << "               Pruned " << num_removed << " non-tight cuts "
	   << "from the LP\n";
      augrounds = 0;
      continue;
    }

    rval = LPcore.pivot_back();
    if(rval) goto CLEANUP;

    cut_rval = CutControl.primal_sep(augrounds, stat);
    if(cut_rval){
      rval = 1;
      goto CLEANUP;
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

    if(cut_rval == 2)
      break;
  }

  if(stat != 3 || rounds == roundlimit){
    cout << "\n Terminated due to: ";
    if(cut_rval == 2){
      cout << "lack of cutting planes.\n      ";
      print.pivot(stat);
      cout << "\n";
    }
    else if(rounds == roundlimit)
      cout << "artificial round limit\n";
  }
  
  cout << "  "
       << (LPcore.numrows() - LPcore.best_tour_nodes.size())
       << " cutting planes added over "
       << rounds << " rounds of separation" << endl;

  CutControl.profile();

    
  cout << "\n Total time for Purecut::solve: "
       << (PSEP_zeit() - routine_start) << "s\n";
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
