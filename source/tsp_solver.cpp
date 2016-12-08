#include<iostream>

#include "tsp_solver.hpp"
#include "lp.hpp"
#include "util.hpp"

using std::vector;
using std::unique_ptr;

namespace CMR {

TSPSolver::TSPSolver(const std::string &fname, RandProb &randprob,
		     OutPrefs _outprefs, LP::Prefs _prefs,
		     unique_ptr<CCdatagroup> &dat,
		     const bool sparse, const int quadnearest) try :
  outprefs(_outprefs), //TODO: for good const form this should be constructed
  //differently maybe, maybe in parse_args by passing pointer
  GraphGroup(fname, outprefs.probname, randprob, dat, sparse, quadnearest,
	     outprefs.dump_xy),
  BestGroup(GraphGroup.m_graph, GraphGroup.delta, dat, outprefs.probname,
	    randprob.seed, outprefs.save_tour, outprefs.save_tour_edges),
  LPGroup(GraphGroup.m_graph, _prefs, BestGroup.perm),
  PureCut(CMR::make_unique<CMR::PureCut>(GraphGroup, BestGroup, LPGroup,
					   SupportGroup, outprefs))
{
}
catch (...) {
  std::cerr << "Failure in TSPSolver constructor\n";
  throw 1;
 }

//ugly temporary workaround, should fix later
static CMR::RandProb dummy;

TSPSolver::TSPSolver(const std::string &fname, const std::string &tourname,
		     CMR::OutPrefs _outprefs, CMR::LP::Prefs _prefs,
		     std::unique_ptr<CCdatagroup> &dat,
		     const bool sparse, const int quadnearest) try :
  outprefs(_outprefs),
  GraphGroup(fname, outprefs.probname, dummy, dat, sparse,
	     quadnearest, outprefs.dump_xy),
  BestGroup(tourname, GraphGroup.m_graph, GraphGroup.delta, dat,
	    outprefs.probname, outprefs.save_tour, outprefs.save_tour_edges),
  LPGroup(GraphGroup.m_graph, _prefs, BestGroup.perm),
  PureCut(CMR::make_unique<CMR::PureCut>(GraphGroup, BestGroup, LPGroup,
					   SupportGroup, outprefs))
{
}
catch (...) {
  std::cerr << "Failure in TSPSolver constructor\n";
  throw 1;
 }  

int TSPSolver::call(SolutionProtocol solmeth, const bool sparse) {
  LP::PivType piv_status;
  int rval = 0;
  
  if (solmeth == SolutionProtocol::PURECUT) {
    PivotPlan plan;
    if (sparse)
      plan = PivotPlan(GraphGroup.m_graph.node_count, PivPresets::SPARSE);      
   
    if (PureCut->solve(plan, piv_status)) {
      std::cerr << "TSPSolver.call(PURECUT) failed\n";
      return 1;
    }
    
    return 0;
    
  } else {
    PivotPlan plan(GraphGroup.m_graph.node_count, PivPresets::ROOT);
    if (PureCut->solve(plan, piv_status)) {
      std::cerr << "TSPSolver.call(ABC) failed to call PureCut at root\n";
      return 1;
    }

    if (piv_status == LP::PivType::FathomedTour) return 0;
    if (piv_status == LP::PivType::Subtour) {
      std::cerr << "Terminated with inseparable subtour inequality\n";
      return 1;
    }

    int ecount = GraphGroup.m_graph.edge_count;
    std::vector<double> lower_bounds;
    try{ lower_bounds.resize(ecount); }
    catch (const std::bad_alloc &) {
      rval = 1; CMR_GOTO_CLEANUP("Couldn't allocate lower bounds, ");
    }

    if (CMRlp_getlb(&(LPGroup.m_lp), &lower_bounds[0], 0, ecount - 1)) {
      std::cerr << "TSPSolver.call(ABC) failed to get lower bounds\n";
      return 1;
    }

    ABC.reset(new CMR::ABC(BB::BranchPlan::Main,
			    GraphGroup, BestGroup, LPGroup, SupportGroup,
			    lower_bounds, *PureCut));

    return ABC->solve();
  }

 CLEANUP:
  if (rval)
    std::cerr << "TSPSolver.call failed\n";
  return rval;
}
  
}
