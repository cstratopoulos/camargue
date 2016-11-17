/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief LIGHT SIMPLE DP SEPARATION
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_SIMPLEDP_HPP
#define PSEP_SIMPLEDP_HPP

#include "cuts.hpp"
#include "tooth.hpp"
#include "Graph.hpp"

#include <vector>

namespace PSEP {

/** Class for light simple DP separation.
 * Given a tour and an lp solution, this class will perform exact primal
 * separation of light simple domino parity (DP) inequalities, using the 
 * standard algorithm of Fleischer, Letchford, and Lodi (2006) with the 
 * primal candidate approach in Letchford-Lodi (Primal Separation Algorithms),
 * and some enhancements from my thesis. 
 */
template<> class Cut<dominoparity> {
public:
  Cut<dominoparity>(std::vector<int> &_delta, std::vector<int> &_edge_marks,
		    std::vector<int> &_best_tour_nodes,
		    std::vector<int> &_perm,
		    SupportGraph &_G_s,
		    std::vector<int> &_support_elist,
		    std::vector<double> &_support_ecap);

  int cutcall();

protected:
  int separate();
  int add_cuts();

private:
  int build_external(const dominoparity &dp_cut);

  CandidateTeeth candidates;
};

}

#endif
