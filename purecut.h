#ifndef PSEP_PURE_CUT_H
#define PSEP_PURE_CUT_H

#include <iostream>
#include <vector>

#include "datagroups.h"
#include "Graph.h"
#include "lp.h"
#include "PSEP_util.h"
#include "cutcall.h"
#include "LPcore.h"
#include "printer.h"
#include "augheuristic.h"

class PSEP_PureCut {
 public:
  PSEP_PureCut(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
	       PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup):
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
    cutcall(GraphGroup, BestGroup, LPGroup, SupportGroup),
      LPcore(LPGroup, GraphGroup, SupportGroup, BestGroup),
      Aug(BestGroup, LPGroup, SupportGroup){}

  int solve(const bool heuristic);
  PSEP_Printer print;
  
 private:
  PSEP_Cutcall cutcall;
  PSEP_LP_Core LPcore;
  PSEP_AugHeuristic Aug;
};

#endif
