#include "branch_util.hpp"

#include <algorithm>

using std::vector;

namespace CMR {
namespace ABC {

static inline bool scorepair_greater(const ScorePair &p, const ScorePair &o)
{
    return p.first > o.first;
}

double var_score(double mult, double v0, double v1)
{
    return  ((v0 < v1) ? (mult * v0 + v1) : (v0 + mult * v1)) / (mult + 1);
}

vector<ScorePair> ranked_cands(const vector<int> &cand_inds,
                               const vector<int> &down_est,
                               const vector<int> &up_est,
                               double mult, int num_return)
{
    vector<ScorePair> result;
    
    for (int i = 0; i < cand_inds.size(); ++i)
        result.push_back(ScorePair(cand_inds[i], var_score(mult, down_est[i],
                                                           up_est[i])));

    std::sort(result.begin(), result.end(), scorepair_greater);

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

}
}
