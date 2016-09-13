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

#include "datagroups.h"
#include "Graph.h"
#include "lp.h"
#include "PSEP_util.h"
#include "cutcall.h"
#include "pivplan.h"
#include "LPcore.h"
#include "LPfixing.h"
#include "printer.h"

namespace PSEP {
  class ABC;
  namespace BB {
    class Visitor;
  }
  
  class PureCut {
  public:
  PureCut(Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
	  Data::LPGroup &LPGroup, Data::SupportGroup &SupportGroup):
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
      CutControl(GraphGroup, BestGroup, LPGroup, SupportGroup),
      LPPrune(GraphGroup, LPGroup),
      LPCore(LPGroup, GraphGroup, SupportGroup, BestGroup, LPPrune),
      LPFix(BestGroup, GraphGroup, LPGroup){}


    int solve(PSEP::PivotPlan &plan, PSEP::LP::PivType &piv_stat);
    PSEP_Printer print;
  
  private:
    friend class PSEP::BB::Visitor;
    friend class PSEP::ABC;
    
    PSEP::CutControl CutControl;
    PSEP::LP::CutPrune LPPrune;
    PSEP::LP::Core LPCore;
    PSEP::LP::EdgeFix LPFix;
  };
}

#endif
