#include "cutmon.hpp"
#include "util.hpp"
#include "err_util.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::vector;


namespace CMR {
namespace Eps = Epsilon;

namespace LP {

CutMonitor::CutMonitor(int cutcount) try
    : pivot_ages(cutcount, 0)
{} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("CutMonitor constructor failed");
}

CutMonitor::CutMonitor(CutMonitor &&CM) noexcept
    : pivot_ages(std::move(CM.pivot_ages)) {}


CutMonitor &CutMonitor::operator=(CutMonitor &&CM) noexcept
{
    pivot_ages = std::move(CM.pivot_ages);
    return *this;
}

/**
 * @param[in] piv_duals dual values for each cut in the LP relaxation at the
 * current non-degenerate pivot.
 * @pre \p piv_duals shall be at least as big as pivot_ages.
 * @post pivot_ages is the same size as \p piv_duals. For every entry of
 * \p piv_duals greater than or equal to CMR::Epsilon::Cut, the corresponding
 * entry of pivot_ages is reset to zero. Else, it is incremented by one.
 */
void CutMonitor::update_pivs(const vector<double> &piv_duals)
{
    if (piv_duals.size() < pivot_ages.size())
        throw logic_error("piv_duals too small in CutMonitor::update_pivs");
    else {
        try { pivot_ages.resize(piv_duals.size(), 0.0);}
        catch (const exception &e) {
            cerr << e.what() << endl;
            throw runtime_error("CutMonitor::update_pivs resize failed");
        }
    }

    for (int i = 0; i < piv_duals.size(); ++i)
        if (piv_duals[i] < Eps::Cut)
            ++pivot_ages[i];
        else
            pivot_ages[i] = 0;
}

}
}
