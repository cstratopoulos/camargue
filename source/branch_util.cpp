#include "branch_util.hpp"

#include <algorithm>

using std::vector;

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
    return  ((v0 < v1) ? (mult * v0 + v1) : (v0 + mult * v1)) / (mult + 1);
}

/**
 * @param[in] cand_inds candidate branching indices
 * @param[in] down_est fix-to-zero estimates for \p cand_inds
 * @param[in] up_est fix-to-one estimates for \p cand_inds
 * @param[in] num_return the max number of variables to be returned
 * @pre `cand_inds.size() == down_est.size() == up_est.size()`
 * @pre each entry of \p down_est and \p up_est is the estimate for the
 * corresponding entry of \p cand_inds.
 * @returns A vector of at most \p num_return ScorePair objects, the 
 * \p num_return highest scored indices from \p cand_inds found by 
 * applying var_score to each index.
 */

vector<ScorePair> ranked_cands(const vector<int> &cand_inds,
                               const vector<double> &down_est,
                               const vector<double> &up_est,
                               double mult, int num_return)    
{
    vector<ScorePair> result;
    
    for (int i = 0; i < cand_inds.size(); ++i)
        result.push_back(ScorePair(cand_inds[i], var_score(mult, down_est[i],
                                                           up_est[i])));

    std::sort(result.begin(), result.end(),
              [](const ScorePair &o, const ScorePair &p)
              { return o.first > p.first; });

    result.resize(std::min((int) result.size(), num_return));
    return result;
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
