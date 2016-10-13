#include "purecut.hpp"

#include<iomanip>

using namespace std;
using namespace PSEP;

int PureCut::solve(PivotPlan &plan, LP::PivType &piv_stat){
  int rval = 0, cut_rval;

  double piv_val;
  int rounds = 0, augrounds = 0;

  bool prune = !plan.is_branch();

  double pivtime, total_pivtime = 0, max_pivtime = 0;
  int num_removed = 0;
  double routine_start, fixing_start, routine_total;
  double fixtime = 0;

  // if(plan.perform_elim()){
  //   fixing_start = zeit();
  //   rval = LPFix.redcost_fixing();
  //   if(rval) goto CLEANUP;
  //   fixing_start = zeit() - fixing_start;
  //   fixtime += fixing_start;
  // }
  
  if(!plan.is_branch())
    cout << "Pivoting until optimality or no more cuts" << endl;

  routine_start = zeit();
  plan.start_timer();
  
  while(plan.condition(augrounds)){
    rounds++;
    augrounds++;

    if(rounds % 50 == 0){
      if(plan.perform_elim()){
	cout << "Calling edge elimination again...\n\n  ";
	fixing_start = zeit();
	rval = LPFix.redcost_fixing();
	if(rval) goto CLEANUP;
	fixing_start = zeit() - fixing_start;
	fixtime += fixing_start;

	plan.current_edge_ratio = LPCore.numcols() / plan.ncount;

	rval = LPCore.rebuild_basis(prune);
	if(rval) goto CLEANUP;
      }
    }

    pivtime = zeit();
    rval = LPCore.pivot_until_change(piv_stat);
    if(rval) goto CLEANUP;
    pivtime = zeit() - pivtime;
    total_pivtime += pivtime;
    if(pivtime > max_pivtime) max_pivtime = pivtime;

    piv_val = LPCore.get_obj_val();

    if(piv_stat == LP::PivType::FathomedTour){
      cout << "\n\n    ROUND " << rounds << " -- ";
      print.pivot(piv_stat);
      cout << "                Pivot objval: ";
      printf("%.6f\n", LPCore.get_obj_val());
      cout << "                * * * * * * * * * *\n";
      break;
    }

    if(piv_stat == LP::PivType::Tour){
      if(!LPCore.test_new_tour() || plan.is_branch())
	continue;
      
      if(plan.is_branch())
	break;

      cout << "\n\n    !!!AUGMENTED TOUR!!!! Round " << rounds << "\n    ";
      print.pivot(piv_stat);
      cout << "                Pivot objval: ";
      printf("%.6f\n", LPCore.get_obj_val());
      
      rval = LPCore.update_best_tour();
      if(rval) goto CLEANUP;
      
      rval = LPPrune.prune_cuts(num_removed);
      if(rval) goto CLEANUP;
      cout << "               Pruned " << num_removed << " non-tight cuts "
	   << "from the LP\n";
      augrounds = 0;
      continue;
    }

    cut_rval = CutControl.primal_sep(augrounds, piv_stat);
    if(cut_rval == 1){
      rval = 1;
      goto CLEANUP;
    }

    //pivot_back used to be here
    
    if(cut_rval == 2){
      if(piv_stat == LP::PivType::Subtour){
	rval = LPCore.add_connect_cut();
	if(rval) goto CLEANUP;
      }

      rval = LPCore.pivot_back();
      if(rval) goto CLEANUP;
      
      cout << "\n    Round " << rounds << ", calling safe GMI sep....";
      rval = CutControl.safe_gomory_sep();
      if(rval == 1) goto CLEANUP;
      if(rval == 2) break;

      if(piv_stat == LP::PivType::Subtour){
	rval = LPCore.del_connect_cut();
	if(rval) goto CLEANUP;
      }
    } else { //TODO: this is a bit ungraceful
      rval = LPCore.pivot_back();
      if(rval) goto CLEANUP;
    }

    rval = CutControl.add_primal_cuts(); //TODO: this should add the gomory cuts
    if(rval) goto CLEANUP;

    if(rounds % 50 == 0 && !plan.is_branch()){
      cout << "\n PIVOTING ROUND: " << rounds << " [ "
	   << (LPCore.numrows() - LPCore.best_tour_nodes.size())
	   << " cuts in the LP ]\n";
      cout << "   ";
      print.pivot(piv_stat);
      cout << "   Pivot objval: " << piv_val << "\n";
      cout << "   Avg piv time: " << setprecision(2)
	   << ((double) (total_pivtime / rounds)) <<"s, (longest "
	   << max_pivtime << "s)\n" << setprecision(6);
      max_pivtime = 0;	   
    }
  }

  routine_total = zeit() - routine_start;

  if(plan.is_branch() && piv_stat == LP::PivType::Tour){
    cout << "   Terminated due to augmented tour in branch solve ("
	 << rounds << " rounds)\n";
  } else if(piv_stat != LP::PivType::FathomedTour){
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
	 << (LPCore.numrows() - LPCore.best_tour_nodes.size())
	 << " cutting planes added over "
	 << rounds << " rounds of separation" << endl;

  if(!plan.is_branch())
    CutControl.profile(routine_total);

  cout << "             Total time (s) pivoting: " << total_pivtime
       << ", ratio: "
       << (total_pivtime / routine_total) << "\n";

    
  cout << "\n Total time (s) for Purecut::solve: "
       << (routine_total) << "\n";
  if(plan.perform_elim())
    cout <<"         LPFix::redcost_fixing: "
	 << fixtime << "s\n";

  if(piv_stat != LP::PivType::FathomedTour && plan.perform_elim()){
    LPFix.redcost_fixing();
    LPCore.set_support_graph();
  }


 CLEANUP:
  if(rval)
    cerr << "Error entry point: PureCut::solve()\n";
  return rval;
}
