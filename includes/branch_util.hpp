/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief UTILITY FUNCTIONS AND ENUMS FOR BRANCHING
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_UTIL_H
#define CMR_BRANCH_UTIL_H

#include <utility>
#include <vector>

namespace CMR {

/** Augment-Branch-Cut solution. */
namespace ABC {

/** Orders of magnitude factor for fixing a variable to one. */
constexpr int one_factor = 5;

/** Orders of magnitude factor for fixing a variable to zero. */
constexpr int zero_factor = 10;

/** Protocols for enforcing branching on variables. */
enum class EnforceStrat {
    /** Tighten bounds in the lp relaxation.
     * This protocol enforces branching on variables by adding the literal
     * equality bound.
     */
    Clamp,

    /** Perturb objective function coefficients.
     * When branching in affirmation with the current tour, this protocol
     * makes the corresponding variable agree with the current tour. Otherwise,
     * it sets the objective function to a value several orders of magnitude
     * larger than the current upper bound, hopefully forcing the variable
     * to be chosen or rejected without rendering the tour infeasible.
     */
    Dive,

    /** Add a single row to the relaxation. */
    Naive,

    /** Add a ton of rows to the relaxation. */
    Rows
};

/** Rank a branching variable in terms of its down and up estimates.
 * For some estimation metric, let \p v0 be the estimate obtained for fixing
 * a variable to zero, and let \p v1 be the estimate for fixing the same 
 * variable to one. This function returns \f[
 * \frac{\gamma\min(v_0, v_1) + \max(v_0, v_1)}{\gamma + 1} \f]
 * where \f$ \gamma \f$ is \p mult.
 */
double var_score(double mult, double v0, double v1);

/** Alias declaration for ranking branching variables. 
 * The first entry represents a candidate branching variable index, the second
 * entry represents a priority ranking for branching on that variable, where
 * larger values are better.
 */
using ScorePair = std::pair<int, double>;

/** Produce a list of ranked scored branching edges.
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
std::vector<ScorePair> ranked_cands(const std::vector<int> &cand_inds,
                                    const std::vector<double> &down_est,
                                    const std::vector<double> &up_est,
                                    double mult, int num_return);

int num_digits(const double val); /**< The number of base 10 digits in val. */

/** Get coefficients for imposing objective function diving.
 * @param[in] upper_bound an upper bound on the minimization problem
 * @param[out] zero_coeff objective function coeff for fixing a var to zero
 * @param[out] one_coeff objective function coeff for fixing a var to one.
 * @returns the pair of coeffs \p zero_coeff and \p one_coeff that can be used
 * to implicitly force or exclude a variable by objective function 
 * perturbation.
 */
void dive_coeffs(const double upper_bound,
                 double &zero_coeff, double &one_coeff);

/** Determine if an objective value implies compliance with Dive branching.
 * @param[in] zero_coeff the coefficient used to fix variables to zero
 * @param[in] one_coeff the coefficient used to fix variables to one
 * @param[in] branch_objval the objective value of a solution in the branch
 * tree.
 * @param[in] zero_count the number of variables being fixed to zero.
 * @param[in] one_count the number of variables being fixed to one. 
 * @returns True iff \p branch_objval implies that variables are being fixed
 * to their desired values, false otherwise. 
 */
bool branch_compliant(const double zero_coeff, const double one_coeff,
                      const double branch_objval,
                      const int zero_count, const int one_count);


}
}

#endif
