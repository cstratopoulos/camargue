/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Functions, structs/enums, and constants for branching.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_UTIL_H
#define CMR_BRANCH_UTIL_H

#include "graph.hpp"
#include "lp_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

namespace CMR {

/// Augment-Branch-Cut solution.
namespace ABC {

/// Constants related to branch selection.
namespace SB {

constexpr double LongMult = 10;
constexpr double StrongMult = 100;

constexpr int Cands1 = 5;
constexpr int Cands2 = 2;

constexpr int Lim1Min = 10;
constexpr int Lim2Max = 500;

inline int round1_limit(int avg_itcount)
{
    return std::min(std::max(Lim1Min, 2 * avg_itcount), Lim2Max);
}

inline int round2_limit(int avg_itcount)
{
    return std::min(4 * avg_itcount, Lim2Max);
}

}

/// Return a "large" edge length relative to the capacities in \p ecap.
int large_len(int ncount, const std::vector<Graph::Edge> &edges);

/// Rank a branching variable in terms of its down and up estimates.
double var_score(double mult, double v0, double v1);


/// Get a list of candidate branch edges using fractional long edge branching.
std::vector<int> length_weighted_cands(const std::vector<Graph::Edge> &edges,
                                       const std::vector<int> &indices,
                                       const std::vector<double> &x,
                                       const int num_return);


struct ScoreTuple {
    ScoreTuple() = default;

    ScoreTuple(EndPts e, LP::Estimate down, LP::Estimate up,
               LP::Basis base, double mult, double ub);

    EndPts ends;

    LP::Estimate down_est;
    LP::Estimate up_est;

    double score;

    LP::Basis contra_base;
};

inline std::ostream &operator<<(std::ostream &os, const ScoreTuple &T)
{
    os << T.ends << "\n\tdown looks " << T.down_est.sol_stat << ", val "
       << T.down_est.value << "\n\t";
    os << "up looks " << T.up_est.sol_stat << ", val " << T.up_est.value;
    os << "\n\tscore " << T.score;
    return os;
}

/// Produce a list of fixed max size containing ranked scored branching edges.
std::vector<ScoreTuple> ranked_cands(const std::vector<int> &cand_inds,
                                     std::vector<LP::Estimate> &down_est,
                                     std::vector<LP::Estimate> &up_est,
                                     const std::vector<Graph::Edge> &edges,
                                     std::vector<LP::Basis> &contra_bases,
                                     double mult, double ub, int num_return);



}
}

#endif
