/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Monitoring cuts in the LP relaxation.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CUTMON_H
#define CMR_CUTMON_H

#include "lp_util.hpp"


#include <vector>

namespace CMR {
namespace LP {

/// After how many rounds is a cut age considered old.
/// A cut will be moved to the pool if its dual value has been below
/// CMR::Epsilon::Cut for this many rounds or more.
constexpr int PivotOld = 25;

/// Monitoring the dual values of cuts at successive pivots.
class CutMonitor {
public:
    CutMonitor() = default;

    /// Construct a CutMonitor for an LP relaxation with \p cutcount cuts.
    CutMonitor(int cutcount);

    CutMonitor(CutMonitor &&CM) noexcept;
    CutMonitor &operator=(CutMonitor &&CM) noexcept;

    /// Update pivot_ages with the dual values in \p piv_duals.
    void update_pivs(const std::vector<double> &piv_duals);

    /// Delete an indicated set of cuts.
    void del_cuts(const std::vector<int> &delset, int ncount);

    const std::vector<int> &get_piv_ages() const { return pivot_ages; }

private:
    /// A vector of integers for how old each cut is.
    /// We maintain that pivot_ages is a vector of length equal to the number
    /// of cuts in the LP, representing the number of rounds since a given cut
    /// received a dual variable above CMR::Epsilon::Cut.
    std::vector<int> pivot_ages;
};

}
}

#endif
