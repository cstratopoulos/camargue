#include "branch_util.hpp"
#include "graph.hpp"
#include "util.hpp"

#include <algorithm>
#include <iostream>

using std::vector;

using std::ostream;
using std::cout;
using std::cerr;

namespace CMR {
namespace ABC {

using Ptype = Problem::Type;

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

ScoreTuple::ScoreTuple(int ind, ScorePair down, ScorePair up, double mult,
                       double ub)
    : index(ind), score_priority(std::max(down.first, up.first)),
      down_est(down), up_est(up),
      score(var_score(mult, down_est.second, up_est.second))
{}

bool operator>(ScoreTuple s, ScoreTuple t)
{
    return
    std::tie(s.score_priority, s.score) > std::tie(t.score_priority, t.score);
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
                                const vector<ScorePair> &down_est,
                                const vector<ScorePair> &up_est,
                                const double mult, const double ub,
                                const int num_return)
{
    vector<ScoreTuple> result;

    for (int i = 0; i < cand_inds.size(); ++i)
        result.emplace_back(cand_inds[i], down_est[i], up_est[i], mult, ub);

    std::sort(result.begin(), result.end(), std::greater<ScoreTuple>());

    if (result.size() <= num_return)
        return result;
    else
        return vector<ScoreTuple>(result.begin(), result.begin() + num_return);

}

int num_digits(const double val)
{
    int factor = 10;
    int result = 1;
    
    int check = std::abs(val);

    while ((check % factor) != val) {
        ++result;
        factor *= 10;
    }

    return result;
}


/**
 * @param[in] upper_bound an upper bound on the minimization problem
 * @param[out] zero_coeff objective function coeff for fixing a var to zero
 * @param[out] one_coeff objective function coeff for fixing a var to one.
 * @returns the pair of coeffs \p zero_coeff and \p one_coeff that can be used
 * to implicitly force or exclude a variable by objective function 
 * perturbation.
 */
void dive_coeffs(const double upper_bound,
                 double &zero_coeff, double &one_coeff)
{
    double result = 1;

    int num_dig = num_digits(upper_bound);
    int one_pow = num_dig + one_factor;
    int zero_pow = num_dig + zero_factor;

    int i = 0;

    while (i++ < one_pow)
        result *= 10;

    one_coeff = -result;

    while (i++ < zero_pow)
        result *= 10;

    zero_coeff = result;
}

ostream &operator<<(ostream &os, Ptype type)
{
    switch (type) {
    case Ptype::Root:
        os << "Root";
        break;
    case Ptype::Affirm:
        os << "Affirm";
        break;
    case Ptype::Contra:
        os << "Contra";
        break;
    }

    return os;
}

ostream &operator<<(ostream &os, const Problem &prob)
{
    os << prob.type << " branch on edge "
       << prob.edge_ind;
    return os;
}

/** 
 * @param[in] zero_coeff the coefficient used to fix variables to zero
 * @param[in] one_coeff the coefficient used to fix variables to one
 * @param[in] branch_objval the objective value of a solution in the branch
 * tree.
 * @param[in] zero_count the number of variables being fixed to zero.
 * @param[in] one_count the number of variables being fixed to one. 
 * @returns True iff \p branch_objval implies that variables are being fixed
 * to their desired values, false otherwise. 
 */

}
}
