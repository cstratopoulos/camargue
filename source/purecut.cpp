#include "purecut.hpp"

#include<iomanip>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setprecision;

namespace PSEP {

int PureCut::solve(PivotPlan &plan, LP::PivType &piv_stat){
  int rval = 0, cut_rval = 0;

  double piv_val;
  int rounds = 0, augrounds = 0;

  bool prune = !plan.is_branch(), haveslack = false;

  Timer pivtimer("Pivot time", &pctime);
  int num_removed = 0;
  double fixing_start;
  double fixtime = 0;
  
  if(!plan.is_branch())
    cout << "Pivoting until optimality or no more cuts" << endl;

  pctime.start();
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

    pivtimer.resume();
    rval = LPCore.pivot_until_change(piv_stat);
    if(rval) goto CLEANUP;
    pivtimer.stop();

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

      haveslack = false;
      
      rval = LPCore.update_best_tour();
      if(rval) goto CLEANUP;
      
      rval = LPPrune.prune_cuts(num_removed);
      if(rval) goto CLEANUP;
      cout << "               Pruned " << num_removed << " non-tight cuts "
	   << "from the LP\n";
      augrounds = 0;
      continue;
    }

    // if(augrounds == 5){
    //   print.lp_edges();
    //   print.best_tour_nodes();
    // }

    cut_rval = CutControl.primal_sep(augrounds, piv_stat);
    if(cut_rval == 1){
      rval = 1;
      goto CLEANUP;
    }

    //pivot_back used to be here
    
    if(cut_rval == 2){
      if(piv_stat == LP::PivType::Subtour){
	LP::PivType connect_stat;
	rval = LPCore.add_connect_cuts(connect_stat);
	if(rval) goto CLEANUP;

	if(connect_stat == LP::PivType::Tour){
	  cout << "Round " << rounds << "\n";
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
      }

      rval = LPCore.pivot_back();
      if(rval) goto CLEANUP;
      
      cout << "\n    Round " << rounds << ", calling safe GMI sep....";
      rval = CutControl.safe_gomory_sep();
      if(rval == 1) goto CLEANUP;
      if(rval == 2){ 
	if(piv_stat == LP::PivType::Subtour){
	  cout << "    Rebuilding basis to accomodate slack subtour\n";
	  rval = LPCore.rebuild_basis(false);
	  if(rval) goto CLEANUP;

	  haveslack = true;

	  continue;
	}
	break;
      }

      if(piv_stat == LP::PivType::Subtour || haveslack){
	rval = LPCore.del_connect_cut();
	if(rval) goto CLEANUP;
	haveslack = false;
      }
    } else { /** @todo this is ungraceful, get rid of the GMI subtour thing. */
      rval = LPCore.pivot_back();
      if(rval) goto CLEANUP;
    }

    rval = CutControl.add_primal_cuts(); /** @todo should add gomory cuts */
    if(rval) goto CLEANUP;

    if(rounds % 50 == 0 && !plan.is_branch()){
      cout << "\n PIVOTING ROUND: " << rounds << " [ "
	   << (LPCore.numrows() - LPCore.best_tour_nodes.size())
	   << " cuts in the LP ]\n";
      cout << "   ";
      print.pivot(piv_stat);
      cout << "   Pivot objval: " << piv_val << "\n";
    }
  }

  pctime.stop();

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
    CutControl.profile();

  pivtimer.report(false);
  pctime.report(true);
  cout << "\n";
  
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

}
