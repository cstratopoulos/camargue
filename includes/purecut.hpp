/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                   PURE CUTTING PLANE SOLVER CLASS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_PURE_CUT_H
#define PSEP_PURE_CUT_H

#include <iostream>
#include <vector>

#include "lp.hpp"
#include "datagroups.hpp"
#include "Graph.hpp"
#include "PSEP_util.hpp"
#include "cutcall.hpp"
#include "pivplan.hpp"
#include "LPcore.hpp"
#include "LPfixing.hpp"
#include "printer.hpp"

namespace PSEP {
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
	  PSEP::OutPrefs &_outprefs):
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
    CutControl(GraphGroup, BestGroup, LPGroup, SupportGroup),
    LPPrune(GraphGroup, LPGroup),
    LPCore(LPGroup, GraphGroup, SupportGroup, BestGroup, LPPrune, _outprefs),
    LPFix(BestGroup, GraphGroup, LPGroup){}


  /** The pure cutting plane control flow loop.
   * This function executes the pure cutting plane solution process, as per
   * the PivotPlan specified in \p plan, with the final outcome of the 
   * function stored in the PivType \p final_piv_stat.
   * @returns 0 if successful, 1 if failure.
   */
  int solve(PSEP::PivotPlan &plan, PSEP::LP::PivType &final_piv_stat);
  
  PSEP_Printer print;
  
private:
  friend class PSEP::BB::Visitor;
  friend class PSEP::ABC;

  PSEP::CutControl CutControl; /**< Manages separation routines/cut adding. */
  PSEP::LP::CutPrune LPPrune; /**< Manages the pruning of cuts from the LP */
  PSEP::LP::Core LPCore; /**< Carries out fundamental LP operations */
  PSEP::LP::EdgeFix LPFix; /**< Fixes edges by reduced cost. */
};
}

#endif
