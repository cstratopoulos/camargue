#include "purecut.h"

#include<iomanip>

using namespace std;
using namespace PSEP;

int PureCut::solve(PivotPlan &plan, PivType &piv_stat){
  int rval = 0, cut_rval;

  double piv_val;
  int rounds = 0, augrounds = 0;

  bool prune = !plan.is_branch();

  double pivtime, total_pivtime = 0, max_pivtime = 0;
  int num_removed = 0;
  double routine_start, fixing_start;

  if(plan.perform_elim()){
    fixing_start = PSEP_zeit();
    rval = LPfix.redcost_fixing();
    if(rval) goto CLEANUP;
    fixing_start = PSEP_zeit() - fixing_start; 
  }
  
  rval = LPcore.basis_init();
  if(rval) goto CLEANUP;

  if(!plan.is_branch())
    cout << "Pivoting until optimality or no more cuts" << endl;

  routine_start = PSEP_zeit();
  plan.start_timer();
  
  while(plan.condition(augrounds)){
    rounds++;
    augrounds++;

    if(rounds % 50 == 0){
      if(plan.perform_elim()){
	cout << "Calling edge elimination again...\n\n  ";
	rval = LPfix.redcost_fixing();
	if(rval) goto CLEANUP;

	plan.current_edge_ratio = LPcore.numcols() / plan.ncount;

	rval = LPcore.rebuild_basis(prune);
	if(rval) goto CLEANUP;
      }
    }

    pivtime = PSEP_zeit();
    rval = LPcore.pivot_until_change(piv_stat);
    if(rval) goto CLEANUP;
    pivtime = PSEP_zeit() - pivtime;
    total_pivtime += pivtime;
    if(pivtime > max_pivtime) max_pivtime = pivtime;

    piv_val = LPcore.get_obj_val();

    if(piv_stat == PivType::FATHOMED_TOUR){
      cout << "\n\n    ROUND " << rounds << " -- ";
      print.pivot(piv_stat);
      cout << "                Pivot objval: "
	   << LPcore.get_obj_val() << "\n";      
      break;
    }

    if(piv_stat == PivType::TOUR){
      if(!LPcore.test_new_tour() || plan.is_branch())
	continue;
      
      if(plan.is_branch())
	break;

      cout << "\n\n    !!!AUGMENTED TOUR!!!!" << endl;
      print.pivot(piv_stat);
      cout << "                Pivot objval: "
	   << LPcore.get_obj_val() << "\n";
      
      rval = LPcore.update_best_tour();
      if(rval) goto CLEANUP;
      
      rval = LPPrune.prune_cuts(num_removed);
      if(rval) goto CLEANUP;
      cout << "               Pruned " << num_removed << " non-tight cuts "
	   << "from the LP\n";
      augrounds = 0;
      continue;
    }

    rval = LPcore.pivot_back();
    if(rval) goto CLEANUP;

    cut_rval = CutControl.primal_sep(augrounds, piv_stat);
    if(cut_rval == 1){
      rval = 1;
      goto CLEANUP;
    }

    if(rounds % 25 == 0 && !plan.is_branch()){
      cout << "\n PIVOTING ROUND: " << rounds << " [ "
	   << (LPcore.numrows() - LPcore.best_tour_nodes.size())
	   << " cuts in the LP ]\n";
      cout << "   ";
      print.pivot(piv_stat);
      cout << "   Pivot objval: " << piv_val << "\n";
      cout << "   Avg piv time: " << setprecision(2)
	   << ((double) (total_pivtime / rounds)) <<"s, (longest "
	   << max_pivtime << "s)\n" << setprecision(6);
      max_pivtime = 0;
	   
    }

    if(cut_rval == 2){
      //      break;
      cout << "\n  Round " << rounds
      	   << ", calling general sep,"
	   << " piv val: " << piv_val << "\n";
      print.pivot(piv_stat);
      cut_rval = CutControl.general_sep(piv_val);
      if(cut_rval == 1) goto CLEANUP;
      if(cut_rval == 2) break;
      cut_rval = 0;
    }
  }

  if(plan.is_branch() && piv_stat == PivType::TOUR){
    cout << "   Terminated due to augmented tour in branch solve ("
	 << rounds << " rounds)\n";
  } else if(piv_stat != PivType::FATHOMED_TOUR){
    cout << "\n Terminated in " << rounds << " rounds due to: ";
    if(!plan.condition(augrounds)){
      plan.profile(augrounds);
    } else if (cut_rval == 2){
      cout << "lack of cutting planes.\n    ";
      print.pivot(piv_stat);
      cout << "\n";
    } else {
      cout << "For uncaught reason?\n";
    }
  }
  
  if(!plan.is_branch())
    cout << "  "
	 << (LPcore.numrows() - LPcore.best_tour_nodes.size())
	 << " cutting planes added over "
	 << rounds << " rounds of separation" << endl;

  if(!plan.is_branch())
    CutControl.profile();

    
  cout << "\n Total time for Purecut::solve: "
       << (PSEP_zeit() - routine_start) << "s\n";
  if(plan.perform_elim())
    cout <<"         LPfix::redcost_fixing: "
	 << fixing_start << "s\n";

  if(piv_stat != PivType::FATHOMED_TOUR && plan.perform_elim()){
    LPfix.redcost_fixing();
    LPcore.set_support_graph();
  }


 CLEANUP:
  if(rval)
    cerr << "Error entry point: PureCut::solve()\n";
  return rval;
}
