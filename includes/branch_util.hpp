/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief UTILITY FUNCTIONS AND ENUMS FOR BRANCHING
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_UTIL_H
#define CMR_BRANCH_UTIL_H

#include "graph.hpp"

#include <utility>
#include <vector>

namespace CMR {

/// Augment-Branch-Cut solution.
namespace ABC {

/// Multiplicity for ranking Driebek penalties.
constexpr double InitialMult = 10;

/// Multiplicity for strong branch ranking.
constexpr double StrongMult = 100;

/// Number of candidates to which we apply 1st round strong branching.
constexpr int SB1Cands = 5;

/// Number of candidates to which we apply 2nd round strong branching.
constexpr int SB2Cands = 2;


/// Iteration limit for first round of strong branching.
constexpr int SB1Lim = 100;

/// Iteration limit for second round of strong branching. 
constexpr int SB2Lim = 500;


constexpr int one_factor = 5; //<! Magnitude factor for fixing var to one.
constexpr int zero_factor = 10; //<! Magnitude factor for fixing var to zero.

/// Rank a branching variable in terms of its down and up estimates.
double var_score(double mult, double v0, double v1);


/// A POD struct for ranking branching edges.
struct ScoreTuple {
    ScoreTuple() = default;

    /// Store the down/up estimates and score for a variable.
    ScoreTuple(int ind, double down, double up, double mult, double ub);
    
    int index; //<! The index of the edge being scored.
    
    double down_est; //<! The estimate for clamping to zero.
    double up_est; //<! The estimate for clamping to one.
    double score; //<! The score composed from down_est and up_est.
};

/// Get a list of candidate branch edges using the J\"unger et al. metric.
std::vector<int> length_weighted_cands(const std::vector<Graph::Edge> &edges,
                                       const std::vector<int> &indices,
                                       const std::vector<double> &x,
                                       const int num_return);

/// Produce a list of fixed max size containing ranked scored branching edges.
std::vector<ScoreTuple> ranked_cands(const std::vector<int> &cand_inds,
                                     const std::vector<double> &down_est,
                                     const std::vector<double> &up_est,
                                     const double mult,
                                     const double ub, const int num_return);

int num_digits(const double val); /**< The number of base 10 digits in val. */

/// Get coefficients for imposing objective function diving.
void dive_coeffs(const double upper_bound,
                 double &zero_coeff, double &one_coeff);

/// Determine if an objective value implies compliance with Dive branching.
bool branch_compliant(const double zero_coeff, const double one_coeff,
                      const double branch_objval,
                      const int zero_count, const int one_count);


}
}

#endif
