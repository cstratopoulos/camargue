/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Managing Core LP relaxations of TSP instances.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CORE_LP_H
#define CMR_CORE_LP_H

#include "lp_interface.hpp"
#include "active_tour.hpp"
#include "cc_lpcuts.hpp"
#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "util.hpp"

#include <algorithm>
#include <utility>
#include <vector>

namespace CMR {

class Solver;

namespace LP {


/** Class for storing the core lp associated to a TSP instance and pivoting.
 * This class contains the edges and constraints currently under consideration
 * in an lp relaxation of a TSP instance.
 */
class CoreLP : public Relaxation {
public:
    /// Construct a CoreLP from combinatorial TSP instance + best tour data.
    CoreLP(Graph::CoreGraph &core_graph_, Data::BestGroup &best_data_);

    /// Compute a primal non-degenerate pivot from the active_tour.
    LP::PivType primal_pivot();

    void pivot_back(bool prune_slacks); //!< Pivot back to active_tour.

    /**@name Cut addition routines.
     * Pop cuts from a queue, adding them to the LP until none remain.
     */
    ///@{

    void add_cuts(Sep::LPcutList &cutq);
    void add_cuts(Sep::CutQueue<Sep::dominoparity> &dp_q);
    void add_cuts(Sep::CutQueue<LP::SparseRow> &gmi_q);
    void add_cuts(Sep::CutQueue<Sep::ex_blossom> &ex2m_q);
    void add_cuts(Sep::CutQueue<Sep::HyperGraph> &pool_q);

    ///@}

    /// Add the edges in \p add_batch to the LP, modifying the core_graph.
    void add_edges(const std::vector<Graph::Edge> &add_batch,
                   bool reinstate);

    /// Remove edges indicated by \p edge_delstat from LP, modifying core_graph.
    void remove_edges(std::vector<int> edge_delstat);


    /// Returns true iff x_vec satisfies all constraints and column bounds.
    template<typename numtype>
    bool check_feas(const std::vector<numtype> &x_vec);

    const Data::SupportGroup &support_data() const { return supp_data; }
    const Sep::ExternalCuts &external_cuts() const { return ext_cuts; }
    const LP::ActiveTour &get_active_tour() const { return active_tour; }

    /// Set active_tour from a list of nodes.
    void set_active_tour(std::vector<int> tour_nodes);

    void tourless_mode(); //!< Set active_tour to tourless mode.

    /// Average number of iterations per primal_pivot.
    int avg_itcount() const { return sum_it_count / num_nd_pivots; }

    double active_tourlen() const { return active_tour.length(); }

    /// The global upper bound for this instance.
    double global_ub() const { return best_data.min_tour_value; }

    bool verbose = false;

    friend class CMR::Solver;

private:
    /// Instate active_tour, checking feasibility if it throws.
    void instate_active();

    /// As above but with a reset instate.
    void reset_instate_active();

    /// Update active_tour via an augmenting pivot.
    void handle_aug_pivot(std::vector<int> tour_nodes, Basis aug_base);

    void prune_slacks(); //!< Prune cuts which are not tight at active_tour.

    void purge_gmi(bool instate); //!< Get rid of any GMI cuts in the LP.

    Graph::CoreGraph &core_graph;
    Data::BestGroup &best_data;
    Data::SupportGroup supp_data;

    Sep::ExternalCuts ext_cuts;

    ActiveTour active_tour;

    std::vector<double> lp_edges;
    std::vector<double> pi_vals;
    std::vector<double> feas_stat;

    int num_nd_pivots = 0;
    int sum_it_count = 0;

    int prev_numrows;

    bool steepest_engaged = false;
};


}
}

#endif
