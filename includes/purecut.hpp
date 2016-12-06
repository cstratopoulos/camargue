/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                   PURE CUTTING PLANE SOLVER CLASS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_PURE_CUT_H
#define CMR_PURE_CUT_H

#include <iostream>
#include <vector>

#include "lp.hpp"
#include "datagroups.hpp"
#include "Graph.hpp"
#include "util.hpp"
#include "cutcall.hpp"
#include "pivplan.hpp"
#include "LPcore.hpp"
#include "LPfixing.hpp"
#include "printer.hpp"
#include "timer.hpp"

namespace CMR {
/* Some forward declarations to allow for friend classes to be declared */
class ABC;
namespace BB {
class Visitor;
}

/** Pure primal cutting plane solution class.
 * This class attempts to solve the TSP by pure primal cutting plane method,
 * i.e., by pivoting and adding primal cuts until certifying an optimal tour
 * or until no primal cuts are found. It can be used on its own or embedded
 * into an ABC solution process.
 */
class PureCut {
public:

  /** The basic constructor.
   * The constructor receives one of each type of Data object and uses
   * them to initialize its members for the solution process. All Data groups
   * are assumed to have been initialized, usually by TSPSolver. 
   */
  PureCut(Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
	  Data::LPGroup &LPGroup, Data::SupportGroup &SupportGroup,
	  CMR::OutPrefs &_outprefs):
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
    pctime("PureCut::solve"),
    CutControl(GraphGroup, BestGroup, LPGroup, SupportGroup, &pctime),
    LPPrune(GraphGroup, LPGroup),
    LPCore(LPGroup, GraphGroup, SupportGroup, BestGroup, LPPrune, _outprefs),
    LPFix(BestGroup, GraphGroup, LPGroup){}


  /** The pure cutting plane control flow loop.
   * This function executes the pure cutting plane solution process, as per
   * the PivotPlan specified in \p plan, with the final outcome of the 
   * function stored in the PivType \p final_piv_stat.
   * @returns 0 if successful, 1 if failure.
   */
  int solve(CMR::PivotPlan &plan, CMR::LP::PivType &final_piv_stat);
  
  CMR_Printer print;
  
private:
  friend class CMR::BB::Visitor;
  friend class CMR::ABC;

  CMR::Timer pctime;

  CMR::CutControl CutControl; /**< Manages separation routines/cut adding. */
  CMR::LP::CutPrune LPPrune; /**< Manages the pruning of cuts from the LP */
  CMR::LP::Core LPCore; /**< Carries out fundamental LP operations */
  CMR::LP::EdgeFix LPFix; /**< Fixes edges by reduced cost. */
};
}

#endif
