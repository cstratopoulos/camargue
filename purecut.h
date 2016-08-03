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
  class PureCut {
  public:
  PureCut(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
	       PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup):
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
      CutControl(GraphGroup, BestGroup, LPGroup, SupportGroup),
      LPPrune(GraphGroup, LPGroup),
      LPcore(LPGroup, GraphGroup, SupportGroup, BestGroup, LPPrune),
      LPfix(BestGroup, GraphGroup, LPGroup){}


    int solve(PSEP::PivotPlan &plan);
    PSEP_Printer print;
  
  private:
    PSEP::CutControl CutControl;
    PSEP::LPPrune LPPrune;
    PSEP_LP_Core LPcore;
    PSEP_LPfix LPfix;
  };
}

#endif
