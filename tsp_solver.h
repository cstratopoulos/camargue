#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <vector>
#include <memory>

#include "datagroups.h"
#include "purecut.h"
#include "ABC.h"
#include "pivplan.h"
#include "PSEP_util.h"


namespace PSEP{
  class TSPSolver {
  public:
    TSPSolver(char *fname, PSEP::LP::Prefs _prefs, CCdatagroup *dat);

    int call(PSEP::SolutionProtocol solmeth);
  
  private:
    Data::GraphGroup GraphGroup;
    Data::BestGroup BestGroup;
    Data::SupportGroup SupportGroup;
    Data::LPGroup LPGroup;
  
    std::unique_ptr<PSEP::PureCut> PureCut;
    std::unique_ptr<PSEP::ABC> ABC;
  };
}

#endif


