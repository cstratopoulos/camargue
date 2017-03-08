#include "branch_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::vector;

using std::ostream;
using std::cout;
using std::cerr;

namespace CMR {
namespace ABC {

/**
 * For some estimation metric, evaluate a priority score for branching on a
 * given variable.
 * @param[in] v0 the estimate obtained for fixing a variable to zero
 * @param[in] v1 the estimate obtained for fixing the same variable to one
 * @param[in] mult the weighting factor.
 * @returns \f[
 * \frac{\gamma\min(v_0, v_1) + \max(v_0, v_1)}{\gamma + 1} \f]
 * where \f$ \gamma \f$ is \p mult.
 */
double var_score(double mult, double v0, double v1)
{
    double denom = mult + 1;
    double num;

    if (v0 < v1)
        num = mult * v0 + v1;
    else
        num = v0 + mult * v1;

    return num / denom;
}

constexpr int Int32Max = 2147483647;

/**
 * The intended use is that the negative return value of this function should
 * be assigned as the edge weight to fix edges in a branching tour.
 * @returns a large edge length for an instance with \p ncount nodes and edge
 * weights \p ecap.
 * @remark this is a re-write of a subroutine from build_sparse_dat in
 * getdata.c from Concorde.
 */
int large_len(int ncount, const vector<Graph::Edge> &edges)
{
    int max_len = std::max_element(edges.begin(), edges.end(),
                                   [](const Graph::Edge &e,
                                      const Graph::Edge &f)
                                   { return e.len < f.len; })->len;
    double v = max_len + 1;
    v *= ncount;

    if (256 * v > Int32Max)
        return Int32Max / 256;
    else
        return (max_len + 1) * ncount;
}

/**
 * @param[in] edges the list of Edges in the core lp.
 * @param[in] indices the indices of fractional basic variables from \p edges
 * @param[in] x the lp solution
 * @param[in] num_return the maximum number of candidate indices returned.
 */
vector<int> length_weighted_cands(const vector<Graph::Edge> &edges,
                                  const vector<int> &indices,
                                  const vector<double> &x,
                                  const int num_return)
{
    double max_under = 0.0;
    double min_over = 1.0;

    for (int i : indices) {
        double val = x[i];
        if (val == 0.5) {
            max_under = 0.5;
            min_over = 0.5;
            break;
        }
        if (val < 0.5) {
            if (val > max_under)
                max_under = val;
        } else if (val < min_over)
            min_over = val;
    }

    double lower_bd = 0.75 * max_under;
    double upper_bd = min_over + (0.25 * (1.0 - min_over));
    vector<int> result;

    for (int i : indices)
        if (x[i] >= lower_bd && x[i] <= upper_bd)
            result.push_back(i);

    std::sort(result.begin(), result.end(),
              [&edges] (int a, int b)
              { return edges[a].len > edges[b].len; });

    if (result.size() <= num_return)
        return result;
    else
        return vector<int>(result.begin(), result.begin() + num_return);
}

ScoreTuple::ScoreTuple(EndPts e, LP::Estimate down, LP::Estimate up,
                       LP::Basis base, double mult, double ub)
    : ends(e), down_est(std::move(down)), up_est(std::move(up)),
      score(var_score(mult, down_est.value, up_est.value)),
      contra_base(std::move(base))
{}

/**
 * @param[in] cand_inds the candidate indices to be ranked.
 * @param[in] down_est estimates for setting corresponding entry of cand_inds
 * to zero.
 * @param[in] up_est like down_est but for setting to one.
 * @param[in] mult the multiplier used to generate the variable score.
 * @param[in] ub current upperbound for the problem. If a variable has down or
 * up estimate bigger than this that is a hint the LP was infeasible or can
 * be pruned by lower bound.
 * @param[in] num_return return at most this many ranked candidates.
 * @pre `cand_inds.size() == down_est.size() == up_est.size()`
 * @returns a vector of ScoreTuple objects sorted in non-increasing order by
 * score.
 */
vector<ScoreTuple> ranked_cands(const vector<int> &cand_inds,
                                vector<LP::Estimate> &down_est,
                                vector<LP::Estimate> &up_est,
                                const vector<Graph::Edge> &edges,
                                vector<LP::Basis> &contra_bases,
                                double mult, double ub, int num_return)
{
    vector<ScoreTuple> result;

    for (int i = 0; i < cand_inds.size(); ++i)
        result.emplace_back(edges[cand_inds[i]], std::move(down_est[i]),
                            std::move(up_est[i]),
                            std::move(contra_bases[i]), mult, ub);

    std::sort(result.begin(), result.end(),
              [](const ScoreTuple &a, const ScoreTuple &b)
              { return a.score > b.score; });

    if (result.size() <= num_return)
        return result;
    else {
        result.resize(num_return);
        return result;
    }
}

}
}
