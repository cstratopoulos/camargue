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
      LPfix(BestGroup, GraphGroup, LPGroup){}


    int solve(PSEP::PivotPlan &plan, PSEP::LP::PivType &piv_stat);
    PSEP_Printer print;
  
  private:
    friend class PSEP::BB::Visitor;
    friend class PSEP::ABC;
    
    PSEP::CutControl CutControl;
    PSEP::LPPrune LPPrune;
    PSEP::LP::Core LPCore;
    PSEP_LPfix LPfix;
  };
}

#endif
