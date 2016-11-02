/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                   PURE CUTTING PLANE SOLVER CLASS
 *
 * This is one of the two possible solution protocol classes that may be 
 * instatiated/executed by the TSPSolver class. It engages in a pure primal
 * cutting plane solution protocol, in which we repeatedly pivot from the
 * current best tour and find primal cutting planes, until separation fails
 * or optimality is proved, or the tour is augmented to a better one.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_PURE_CUT_H
#define PSEP_PURE_CUT_H

#include <iostream>
#include <vector>

#include "lp.hpp"
#include "datagroups.hpp"
#include "Graph.hpp"
#include "PSEP_util.hpp"
#include "cutcall.hpp"
#include "pivplan.hpp"
#include "LPcore.hpp"
#include "LPfixing.hpp"
#include "printer.hpp"

namespace PSEP {
  /* Some forward declarations to allow for friend classes to be declared */
  class ABC;
  namespace BB {
    class Visitor;
  }
  
  class PureCut {
  public:
    /*
     * The constructor receives one of every single type of DataGroup
     * These are needed to initiate all its private member classes
     */
  PureCut(Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
	  Data::LPGroup &LPGroup, Data::SupportGroup &SupportGroup,
	  PSEP::OutPrefs &_outprefs):
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
    CutControl(GraphGroup, BestGroup, LPGroup, SupportGroup),
    LPPrune(GraphGroup, LPGroup),
    LPCore(LPGroup, GraphGroup, SupportGroup, BestGroup, LPPrune, _outprefs),
    LPFix(BestGroup, GraphGroup, LPGroup){}


    /*
     * solve is the function for the main pure cutting plane solution loop/
     * control flow. 
     * plan - controls termination condition and cut adding behavior. See
     *    pivplan.h for more info
     * final_piv_stat - the last pivot status upon termination of the solution
     *    loop
     */
    int solve(PSEP::PivotPlan &plan, PSEP::LP::PivType &final_piv_stat);
    PSEP_Printer print;
  
  private:
    friend class PSEP::BB::Visitor;
    friend class PSEP::ABC;

    /*
     * CutControl - manages separation routines and addition of cuts to the LP
     *    see cutcall.h
     * LPPrune - manages the pruning of cuts from the LP, see LPprune.h
     * LPCore - manages core LP functionality, see LPcore.h
     * LPFix - manages fixing/eliminating edges by reduced cost, see LPfixing.h
     */
    PSEP::CutControl CutControl;
    PSEP::LP::CutPrune LPPrune;
    PSEP::LP::Core LPCore;
    PSEP::LP::EdgeFix LPFix;
  };
}

#endif
