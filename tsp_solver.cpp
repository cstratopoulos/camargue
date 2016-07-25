#include "tsp_solver.h"

using namespace std;

TSP_Solver::TSP_Solver(char *fname, PSEP_LP_Prefs _prefs,
		       CCdatagroup *dat) :
  GraphGroup(fname, dat),
  BestGroup(GraphGroup.m_graph, dat),
  LPGroup(GraphGroup.m_graph, _prefs),
  PureCut(GraphGroup, BestGroup, LPGroup, SupportGroup) {}
  

int TSP_Solver::pure_cut(){
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

    if(stat == PIVOT::FATHOMED_TOUR)
      break;

    if(stat == PIVOT::TOUR){
      cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "!!!AUGMENTED TOUR!!!!" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
      cout << "~Call to delete slack cuts should go here~" << endl;
      if(LPcore.update_best_tour())
	goto CLEANUP;
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

 CLEANUP:
  if(rval)
    cerr << "Error entry point: TSP_Solver::pure_cut()\n";
  return rval;
}

int TSP_Solver::simple_test(){
  if(LPcore.basis_init())
    return 1;

  int stat, num_dp = 0, num_bad = 0;
  int num_seg = 0, segval = 0;
  int num_2match = 0, matchval = 0;
  double segtime, matchtime, routine_start;
  int rounds = 0;
  int total_cuts = 0;

  int max_per_round = LPcore.prefs.max_cuts_round;

  bool in_sep = false;

  cout << "Pivoting until solution in subtour polytope...\n";

  while(!in_sep){
    if(LPcore.pivot_until_change(&stat))
     return 1;

    print.pivot(stat);

    if(stat == PIVOT::FATHOMED_TOUR){
      cout << "Somehow solved problem, exiting\n";
      return 0;
    }

    if(stat == PIVOT::TOUR){
      if(LPcore.update_best_tour())
	return 1;
      else
	continue;
    }

    if(LPcore.pivot_back())
      return 1;

    segtime = PSEP_zeit();
    segval = cutcall.segment(max_per_round, &num_seg);
    if(segval == 1)
      return 1;
    segtime = PSEP_zeit() - segtime;

    matchtime = PSEP_zeit();
    matchval = cutcall.blossom(max_per_round - num_seg, &num_2match);
    if(matchval == 1)
      return 1;
    matchtime = PSEP_zeit() - matchtime;

   
    cout << "Added " << num_seg << " segments"
	 << " (in " << segtime << "s)"
	 <<", " << num_2match
	 << " blossoms"
	 << " (in " << matchtime << "s)" << endl;

    total_cuts += num_seg + num_2match;

    if(num_seg == 0 && stat != PIVOT::SUBTOUR){
      if(cutcall.in_subtour_poly(&in_sep))
	return 1;
    }

    if(segval + matchval == 4 && !in_sep){
      cout << "No more cuts to add and still not in subtour polytope :'(\n";
      return 1;
    }
    rounds++;
  }

  cout << "Added " << total_cuts << " cuts in "
       << rounds << " rounds\n";

  // print.best_tour_nodes();
  // print.lp_edges();

  if(in_sep){
    cout << "Solution is in subtour polytope, building collection...\n";
    routine_start = PSEP_zeit();
    cutcall.simpleDP(250, &num_dp, &num_bad);
    cout << num_dp << " dp inequalities added\n";
    cout << num_bad << " bad inequalities found\n";
    cout << (PSEP_zeit() - routine_start) << "s finding candidate teeth "
	 << "and building light cutgraph/GH tree\n";
  }

  return 0;
  
}
