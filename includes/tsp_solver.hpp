/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                MAIN TSP SOLVER CLASS DEFINITION
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TSP_SOLVER_H
#define TSP_SOLVER_H

#include <vector>
#include <memory>

#include "datagroups.hpp"
#include "purecut.hpp"
#include "ABC.hpp"
#include "pivplan.hpp"
#include "PSEP_util.hpp"


namespace PSEP{

/** The master class for controlling and invoking solution protocols.
 * TSPSolver is the controller object for this project. It is responsible 
 * for initializing Data objects and then using them to construct a solution
 * protocol class, either Purecut or ABC. 
 */
class TSPSolver {
public:
  /**  TSPSolver constructor using filename or random problem.
   * ONE OF \p fname or \p randprob may be used to initialize the problem.
   * `if (fname.empty())` then \p randprob will be used, and
   * `if randprob.nodecount == 0` then \p fname will be used. The command line
   *  argument parser is responsible for giving precisely one valid argument.
   * @param[in] fname a TSP instance with `.tsp` suffix using TSPLIB format
   * @param[in] randprob: parameters for generating a random problem; see
   *    PSEP_util.h for information
   * @param[in] _outprefs a PSEP::OutPrefs structure, governing the amount
   * of information to be written to file. 
   * @param[in] prefs a PSEP::LP::Prefs structure, for parameters related to
   * the LP solver and the addition of cuts.
   * @param[in] dat pointer to an uninitialized `CCdatagroup` object
   * @todo put \p sparse and \p quadnearest in an initial prefs struct
   */
  TSPSolver(const std::string &fname, PSEP::RandProb &randprob,
	    PSEP::OutPrefs _outprefs, PSEP::LP::Prefs _prefs,
	    std::unique_ptr<CCdatagroup> &dat,
	    const bool sparse, const int quadnearest);
  
  /** TSPSolver constructor with starting tour from file
   * This constructor works just like the other one, but with an initial tour
   * given by the file \p tourname. There is no support for initializing a tour
   * from file and using a random problem. 
   */
  TSPSolver(const std::string &fname, const std::string &tourname,
	    PSEP::OutPrefs _out_prefs, PSEP::LP::Prefs _prefs,
	    std::unique_ptr<CCdatagroup> &dat,
	    const bool sparse, const int quadnearest);

  /** Invoke the solution process.
   * This function invokes the solution process using the SolutionProtocol
   * specified by \p solmeth.
   * @returns 0 if success, 1 if failure
   */
  int call(PSEP::SolutionProtocol solmeth, const bool sparse);
  
private:
  PSEP::OutPrefs outprefs;
  
  Data::GraphGroup GraphGroup;
  Data::BestGroup BestGroup;
  Data::SupportGroup SupportGroup;
  Data::LPGroup LPGroup;

  std::unique_ptr<PSEP::PureCut> PureCut;
  std::unique_ptr<PSEP::ABC> ABC;
};
}

#endif


