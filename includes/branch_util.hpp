/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Miscellaneous functions, structs/enums, and constants for branching.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_UTIL_H
#define CMR_BRANCH_UTIL_H

#include "lp_util.hpp"

#include <iostream>
#include <utility>
#include <vector>

namespace CMR {

namespace Graph { struct Edge; }

/// Augment-Branch-Cut solution.
namespace ABC {


constexpr double InitialMult = 10; //!< Scale factor for Driebeek penalties.

constexpr double StrongMult = 100; //!< Scale factor for strong branch scores.

constexpr int SB1Cands = 5; //!< Number of 1st round strong branch candidates.

constexpr int SB2Cands = 2; //!< Number of 2nd round strong branch candidates.

constexpr int SB1Lim = 100; //!< 1st round strong branch iteration limit.

constexpr int SB2Lim = 500; //!< 2nd round strong branch iteration limit.


constexpr int one_factor = 5; //!< Magnitude factor for fixing var to one.
constexpr int zero_factor = 10; //!< Magnitude factor for fixing var to zero.

/// Alias declaration for integer ranking and objective value estimate.
/// Higher is better for both entries, to be sorted lexicographically.
using ScorePair = std::pair<int, double>;

/// A simple structure for recording the status of branching subproblems.
struct Problem {
    enum class Type {
        Root, //!< The root lp.
        Affirm, //!< Enforcing agreement with current tour.
        Contra //!< Enforcing departure from current tour.
    };

    Problem() = default;

    Problem(int ind, ScorePair r);

    Problem(int ind, ScorePair r, LP::Basis cbase);

    int edge_ind;
    Type type;
    ScorePair rank;
    LP::Basis::Ptr contra_base;
};

std::ostream &operator<<(std::ostream &os, Problem::Type type);
std::ostream &operator<<(std::ostream &os, const Problem &prob);

/// Strategies for enforcing Problem::Type::Contra branches.
enum class ContraStrat {
    Fix, /// Change the bounds on an edge in the Relaxation.
    Naive, /// Add a single branch constraint to the Relaxation.
};


/// A POD struct for ranking branching edges.
struct ScoreTuple {
    ScoreTuple() = default;

    ScoreTuple(int ind, ScorePair down, ScorePair up, LP::Basis &&cbase,
               double mult, double ub);

    int index; //!< The index of the edge being scored.

    /// How valuable is the estimate obtained. Higher is better.
    int score_priority;

    ScorePair down_est; //!< The estimate for clamping to zero.
    ScorePair up_est; //!< The estimate for clamping to one.

    LP::Basis contra_base;

    double score; //!< The priority score formed from down_est and up_est.
};

bool operator>(const ScoreTuple &s, const ScoreTuple &t);

/// Rank a branching variable in terms of its down and up estimates.
double var_score(double mult, double v0, double v1);

/// Get a list of candidate branch edges using the J\"unger et al. metric.
std::vector<int> length_weighted_cands(const std::vector<Graph::Edge> &edges,
                                       const std::vector<int> &indices,
                                       const std::vector<double> &x,
                                       const int num_return);

/// Produce a list of fixed max size containing ranked scored branching edges.
std::vector<ScoreTuple> ranked_cands(const std::vector<int> &cand_inds,
                                     const std::vector<ScorePair> &down_est,
                                     const std::vector<ScorePair> &up_est,
                                     std::vector<LP::Basis> &contra_bases,
                                     const double mult,
                                     const double ub, const int num_return);

}
}

#endif
