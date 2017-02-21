/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Monitoring the active tour in the solution process.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_ACTIVE_TOUR_H
#define CMR_ACTIVE_TOUR_H

#include "datagroups.hpp"
#include "lp_interface.hpp"
#include "lp_util.hpp"

#include <utility>
#include <vector>

namespace CMR {
namespace LP {

/** Information about the active tour in a CoreLP.
* This class stores the tour which will be used to generate cuts and perform
* pivots. It should be used in a CoreLP in cooperation with a BestGroup,
* although they need not store the same tour.
* @remark In this documentation and in function names, to "instate" an
* ActiveTour in an LP::Relaxation means that the tour is made to be the
* resident solution in the Relaxation: calls to Relaxation::lp_vec or
* Relaxation::get_objval shall coincide (within Epsilon::Zero) of the values
* stored in the tour.
*/
class ActiveTour {
public:
    ActiveTour() = default;

    /// An ActiveTour for a new CoreLP using the Padberg-Hong approach.
    ActiveTour(const Graph::CoreGraph &graph,
               const Data::BestGroup &best_data);

    /// Construct from an augmenting pivot/dfs island.
    ActiveTour(std::vector<int> tour_nodes_,
               std::vector<double> lp_edges, Basis base,
               double lp_objval,
               const Graph::CoreGraph &graph);

    /// Construct a new ActiveTour from scratch, instating it in \p relax.
    ActiveTour(std::vector<int> tour_nodes_,
               LP::Relaxation &relax,
               const Graph::CoreGraph &graph);

    ActiveTour(ActiveTour &&T) noexcept;
    ActiveTour &operator=(ActiveTour &&T) noexcept;

    void set_basis(Basis new_base) { tour_base = std::move(new_base); }

    /// Instate this tour in \p relax.
    void instate(Relaxation &relax);

    /// Get a new basis for this tour from \p relax, instating it.
    void reset_instate(Relaxation &relax);

    double length() const { return tour_len; }

    const std::vector<int> &nodes() const { return tour_nodes; }
    const std::vector<int> &tour_perm() const { return perm; }

    const std::vector<double> &edges() const { return tour_edges; }

    const Basis &base() const { return tour_base; }

private:
    double tour_len;

    std::vector<int> tour_nodes;
    std::vector<int> perm;

    std::vector<double> tour_edges;

    Basis tour_base;
};

}
}

#endif
