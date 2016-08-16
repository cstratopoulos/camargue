#include<iostream>

#include "tsp_solver.h"
#include "lp.h"

using namespace std;
using namespace PSEP;

TSPSolver::TSPSolver(const string &fname, RandProb &randprob, LP::Prefs _prefs,
		     unique_ptr<CCdatagroup> &dat) :
  GraphGroup(fname, randprob, dat),
  BestGroup(GraphGroup.m_graph, dat),
  LPGroup(GraphGroup.m_graph, _prefs, BestGroup.perm){
  
  PureCut.reset(new PSEP::PureCut(GraphGroup, BestGroup, LPGroup,
				  SupportGroup));
}

int TSPSolver::call(SolutionProtocol solmeth){
  LP::PivType piv_status;
  
  if(solmeth == SolutionProtocol::PURECUT){
    PivotPlan plan;
   
    if(PureCut->solve(plan, piv_status)){
      cerr << "TSPSolver.call(PURECUT) failed\n";
      return 1;
    }
    
    return 0;
    
  } else {
    PivotPlan plan(GraphGroup.m_graph.node_count, PivPresets::ROOT);
    if(PureCut->solve(plan, piv_status)){
      cerr << "TSPSolver.call(ABC) failed to call PureCut at root\n";
      return 1;
    }

    if(piv_status == LP::PivType::FathomedTour) return 0;
    if(piv_status == LP::PivType::Subtour){
      cerr << "Terminated with inseparable subtour inequality\n";
      return 1;
    }

    int ecount = GraphGroup.m_graph.edge_count;
    std::vector<double> lower_bounds(ecount);

    if(PSEPlp_getlb(&(LPGroup.m_lp), &lower_bounds[0], 0, ecount - 1)){
      cerr << "TSPSolver.call(ABC) failed to get lower bounds\n";
      return 1;
    }

    ABC.reset(new PSEP::ABC(GraphGroup, BestGroup, LPGroup, SupportGroup,
			    lower_bounds, *PureCut));

    return ABC->solve();
  }
}
  
