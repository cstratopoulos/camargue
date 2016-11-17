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
#include "datagroups.hpp"

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
		    std::vector<double> &_support_ecap,
		    PSEP::CutQueue<PSEP::dominoparity> &_dp_q)  :
    dp_q(_dp_q),
    candidates(_delta, _edge_marks,
	       _best_tour_nodes, _perm,
	       _G_s, _support_elist, _support_ecap),
    best_tour_nodes(_best_tour_nodes) {}

  int cutcall();
  int parse_cut(const PSEP::dominoparity &dp_cut,
		PSEP::Data::GraphGroup &g_dat,
		PSEP::SupportGraph &G_s, std::vector<double> &rmatval,
		double &rhs);

protected:  
  int separate();
  int add_cuts();

private:

  /** Aggregates the coefficients of a simple tooth inequality.
   * If \p T is the tooth \f$ (i, S) \f$, this function will aggregate the
   * simple tooth inequality 
   * \f[ 2x(E(S)) + x(E(i:S)) \le 2|S| - 1. \f]
   * by adding the lefthand side coefficients to \p rmatval, and the 
   * righthand side to \p rhs.
   */
  void parse_tooth(const PSEP::SimpleTooth *T, PSEP::Data::GraphGroup &g_dat,
		   PSEP::SupportGraph &G_s,
		   std::vector<double> &rmatval, double &rhs);

  /** Aggregates the degree equations for handle nodes.
   * If handle nodes is some subset \f$ H \f$ of vertices, this function
   * will sum the degree equations 
   * \f[ 2x(E(H)) + x(\delta(H)) \le 2|H| \f]
   * into \p rmatval and \p rhs.
   */
  void parse_handle(const std::vector<int> &handle_nodes,
		    PSEP::Data::GraphGroup &g_dat,
		    std::vector<double> &rmatval, double &rhs);

  /** Aggregates the nonnegativity equations for an edge.
   * If \p edge_ends is \f$ (i, j) \f$, this function will aggregate the 
   * nonnegativity inequality \f$ x_{ij} \le 0 \f$ by subtracting one 
   * from entry \f$ (i,j) \f$ of \p rmatval.
   */
  void parse_nonneg_edges(const IntPair &edge_ends,
			  PSEP::Data::GraphGroup &g_dat,
			  std::vector<double> &rmatval);

  PSEP::CutQueue<dominoparity> &dp_q;
  PSEP::CandidateTeeth candidates;
  std::vector<int> &best_tour_nodes;
};

}

#endif
