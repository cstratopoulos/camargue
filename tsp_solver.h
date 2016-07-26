#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <vector>

#include "datagroups.h"
#include "purecut.h"
#include "ABC.h"
#include "PSEP_util.h"



class TSP_Solver {
 public:
  TSP_Solver(char *fname, PSEP_LP_Prefs _prefs, CCdatagroup *dat);

  PSEP_GraphGroup GraphGroup;
  PSEP_BestGroup BestGroup;
  PSEP_SupportGroup SupportGroup;
  PSEP_LPGroup LPGroup;
  
  PSEP_PureCut PureCut;
  PSEP_ABC ABC;

  int call(const bool heur){
    return PureCut.solve(heur);
  }
};

#endif


