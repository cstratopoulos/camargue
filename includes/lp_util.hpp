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

/// Alias declaration for branching estimates.
/// First element: is the problem infeasible.
/// Second element: the objective value estimate.
using InfeasObj = std::pair<bool, double>;

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
