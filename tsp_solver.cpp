#include<iostream>

#include "tsp_solver.h"
#include "lp.h"

using namespace std;

TSP_Solver::TSP_Solver(char *fname, PSEP::LP::Prefs _prefs,
		       CCdatagroup *dat) :
  GraphGroup(fname, dat),
  BestGroup(GraphGroup.m_graph, dat),
  LPGroup(GraphGroup.m_graph, _prefs, BestGroup.perm){
  
  PureCut.reset(new PSEP::PureCut(GraphGroup, BestGroup, LPGroup,
				  SupportGroup));
}

int TSP_Solver::call(PSEP::SolutionProtocol solmeth){
  PSEP::LP::PivType piv_status;
  
  if(solmeth == PSEP::SolutionProtocol::PURECUT){
    PSEP::PivotPlan plan;
   
    if(PureCut->solve(plan, piv_status)){
      cerr << "TSP_Solver.call(PURECUT) failed\n";
      return 1;
    }
    
    return 0;
    
  } else {
    PSEP::PivotPlan plan(GraphGroup.m_graph.node_count, PSEP::PivPresets::ROOT);
    if(PureCut->solve(plan, piv_status)){
      cerr << "TSP_Solver.call(ABC) failed to call PureCut at root\n";
      return 1;
    }

    if(piv_status == PSEP::LP::PivType::FathomedTour) return 0;
    if(piv_status == PSEP::LP::PivType::Subtour){
      cerr << "Terminated with inseparable subtour inequality\n";
      return 1;
    }

    int ecount = GraphGroup.m_graph.edge_count;
    std::vector<double> lower_bounds(ecount);

    if(PSEPlp_getlb(&(LPGroup.m_lp), &lower_bounds[0], 0, ecount - 1)){
      cerr << "TSP_Solver.call(ABC) failed to get lower bounds\n";
      return 1;
    }

    ABC.reset(new PSEP::ABC(GraphGroup, BestGroup, LPGroup, SupportGroup,
			    lower_bounds, *PureCut));

    return ABC->solve();
  }
}
  
