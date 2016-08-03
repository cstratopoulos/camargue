#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <vector>
#include <memory>

#include "datagroups.h"
#include "purecut.h"
#include "ABC.h"
#include "pivplan.h"
#include "PSEP_util.h"



class TSP_Solver {
 public:
  TSP_Solver(char *fname, PSEP_LP_Prefs _prefs, CCdatagroup *dat);

  PSEP_GraphGroup GraphGroup;
  PSEP_BestGroup BestGroup;
  PSEP_SupportGroup SupportGroup;
  PSEP_LPGroup LPGroup;
  
  std::unique_ptr<PSEP::PureCut> PureCut;
  std::unique_ptr<PSEP::ABC> ABC;

  int call(){
    PSEP::PivotPlan plan(GraphGroup.m_graph.node_count,
			 PSEP::PivPresets::ROOT);
    return PureCut->solve(plan);
  }
};

#endif


