#include "tsp_solver.h"

using namespace std;

TSP_Solver::TSP_Solver(char *fname, PSEP_LP_Prefs _prefs,
		       CCdatagroup *dat) :
  GraphGroup(fname, dat),
  BestGroup(GraphGroup.m_graph, dat),
  LPGroup(GraphGroup.m_graph, _prefs, BestGroup.perm)// ,
  // PureCut((PSEP_PureCut *) NULL), ABC((PSEP_ABC *) NULL)
{
  PureCut.reset(new PSEP_PureCut(GraphGroup, BestGroup, LPGroup, SupportGroup));
   //   ABC = new PSEP_ABC(GraphGroup, BestGroup, LPGroup, SupportGroup);
}
  
/*
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
*/
