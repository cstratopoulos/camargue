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

constexpr int PivotOld = 10;

class CutMonitor {
public:
    CutMonitor(int cutcount);
    CutMonitor(CutMonitor &&CM) noexcept;
    CutMonitor &operator=(CutMonitor &&CM) noexcept;

    void update_pivs(const std::vector<double> &piv_duals);

    const std::vector<int> &get_piv_ages() const { return pivot_ages; }
    const std::vector<int> &get_tour_ages() const { return tour_ages; }

private:
    std::vector<int> pivot_ages;
    std::vector<int> tour_ages;

};

}
}

#endif
