/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                MAIN TSP SOLVER CLASS DEFINITION
 *
 * This is the overarching class responsible for managing problem data and 
 * executing different solution protocols. 
 *
 * public:
 *  TSPSolver(string &fname, RandProb &randprob, LP::Prefs prefs,
 *            unique_ptr<CCdatagroup> &dat)
 *  The default constructor for the TSPSolver class.
 *  ONE OF fname or randprob may be used to initialize the problem
 *  fname.empty() implies use of randprob, randprob.nodecount of zero implies
 *  use of fname. The command line argument parser is responsible for 
 *  giving precisely one valid argument. (PSEP.cpp in this case)
 *               -fname: a TSP instance with .tsp suffix using TSPLIB format
 *               -ranndprob: parameters for generating a random problem; see
 *                    PSEP_util.h for information
 *               -prefs: LP solution preferences, see PSEP_util.h for info
 *               -dat: a unique pointer to an uninitialized CCdatagroup object,
 *                    used by Concorde to initialize the problem and also for
 *                    some other Concorde routines. Direct manipulation of  this
 *                    structure is handled entirely by Concorde
 *
 * call(SolutionProtocol solmeth)
 * The function to invoke the TSP solver with solution protocol specified by
 * solmeth. See PSEP_util.h and below for info
 *
 * private:
 *                See datagroups.h for information on all structures in the
 *                Data:: namespace.
 *
 *                PureCut, ABC are pointers to solver protocol objects,
 *                    denoting  either pure cutting plane solution, or an 
 *                    Augment-Branch-Cut approach. An ABC object is initialized
 *                    only if call(solmeth) is called with solmeth ABC, but
 *                    PureCut is initialized with either PURECUT or ABC, as the
 *                    PureCut solver is embedded in the ABC process. 
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
    TSPSolver(const std::string &fname, PSEP::RandProb &randprob,
	      PSEP::LP::Prefs _prefs,
	      std::unique_ptr<CCdatagroup> &dat);

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


