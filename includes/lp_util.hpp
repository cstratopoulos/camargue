/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Miscellaneous functions, structs/enums, and constants for LPs.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_LP_UTIL_H
#define CMR_LP_UTIL_H

#include <iostream>
#include <memory>
#include <utility>
#include <vector>


namespace CMR {

/// Manners related to LP relaxations and solver interfaces.
namespace LP {

/// Constants related to ages of cuts.
/// The age of a cut (wrt a tour or pivot) is for how many dual solutions the
/// cut has assumed a dual value less than Epsilon::DualDust.
namespace CutAge {

constexpr int Babby = -1; //!< A new or reset cut.
constexpr int PivOld = 30; //!< Old cut age at pivot duals.
constexpr int TourOld = 20; //!< Old cut age at tour dual solution.

}

/// Enum class for categorizing lp solutions.
enum class PivType {
    Frac, //!< Fractional solution.
    Subtour, //!< Integral subtour.
    Tour, //!< A new or augmented tour.
    FathomedTour //!< A tour with a dual feasible basis in the current lp.
};

inline bool is_tour_piv(PivType P)
{
    return P == PivType::Tour || P == PivType::FathomedTour;
}

enum BStat {
    AtLower = 0,
    Basic = 1,
    AtUpper = 2,
    FreeSuper = 3
};

struct Basis {
    Basis() = default;
    Basis(Basis &&B) noexcept
        : colstat(std::move(B.colstat)), rowstat(std::move(B.rowstat)) {}

    Basis &operator=(Basis &&B) noexcept
        {
            colstat = std::move(B.colstat);
            rowstat = std::move(B.rowstat);
            return *this;
        }

    bool empty() const { return colstat.empty() || rowstat.empty(); }

    std::vector<int> colstat;
    std::vector<int> rowstat;

    using Ptr = std::unique_ptr<Basis>;
};

/// Struct for storing info from branching estimates.
struct Estimate {
    Estimate() = default;

    Estimate(double val)
        : sol_stat(Stat::Abort), value(val), sb_base()  {}

    Estimate(Estimate &&E) noexcept
        : sol_stat(E.sol_stat), value(E.value), sb_base(std::move(E.sb_base))
        {}

    Estimate &operator=(Estimate &&E) noexcept
        {
            sol_stat = E.sol_stat;
            value = E.value;
            sb_base = std::move(E.sb_base);

            return *this;
        }

    /// Solution statuses from strong branch pivoting.
    enum class Stat {
        Abort, //!< Reached iteration limit.
        Infeas, //!< LP proved infeasible.
        Prune //!< LP optimal with objval that implies pruning.
    } sol_stat;

    double value; //!< The objective value estimate.

    /// The basis from the strong branch search.
    /// If sol_stat is Abort, this is null or the first primal feasible basis
    /// encountered. Otherwise, it is the basis that certifies either of the
    /// other two Stat values.
    Basis::Ptr sb_base;
};

inline std::ostream &operator<<(std::ostream &os, Estimate::Stat estat)
{
    using Estat = Estimate::Stat;
    if (estat == Estat::Prune)
        os << "Prunable";
    else if (estat == Estat::Infeas)
        os << "Infeasible";
    else
        os << "Limit Abort";

    return os;
}

/// Simple struct representing sparse matrix row for passing to LP solver.
struct SparseRow {
    std::vector<int> rmatind; //!< Indices of nonzero entries.
    std::vector<double> rmatval; //!< Coefficients for indices in rmatind.
    char sense = '\0'; //!< 'G' for >=, 'L' for <=, 'E' for ==.
    double rhs = 0.0; //!< The righthand side.
    double lp_viol = 0.0; //!< (Optional) violation wrt some vector.
};

/// Operator overload for writing LP::PivType to output stream.
inline std::ostream &operator<<(std::ostream &os, PivType piv)
{
    using Ptype = LP::PivType;

    if (piv == Ptype::Frac)
        os << "Fractional";
    else if (piv == Ptype::Subtour)
        os << "Integral subtour";
    else if (piv == Ptype::Tour)
        os << "Tour";
    else if (piv == Ptype::FathomedTour)
        os << "Optimal tour";

    return os;
}


}
}

#endif
