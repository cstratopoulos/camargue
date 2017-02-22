/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Managing Core LP relaxations of TSP instances.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CORE_LP_H
#define CMR_CORE_LP_H

#include "lp_interface.hpp"
#include "active_tour.hpp"
#include "cc_lpcuts.hpp"
#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "cutmon.hpp"
#include "util.hpp"

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
    CoreLP(Graph::CoreGraph &core_graph_,
           Data::BestGroup &best_data);

    /// Compute a primal non-degenerate pivot from the active_tour.
    LP::PivType primal_pivot();

    void pivot_back(bool prune_slacks); //!< Pivot back to active_tour.

    void add_cuts(Sep::LPcutList &cutq);
    void add_cuts(Sep::CutQueue<Sep::dominoparity> &dp_q);
    void add_cuts(Sep::CutQueue<LP::SparseRow> &gmi_q);
    void add_cuts(Sep::CutQueue<Sep::ex_blossom> &ex2m_q);
    void add_cuts(Sep::CutQueue<Sep::HyperGraph> &pool_q);

    void add_edges(const std::vector<Graph::Edge> &add_batch);

    /// Get a const reference to the SupportGroup for the most recent pivot.
    const Data::SupportGroup &support_data() const { return supp_data; }

    const Sep::ExternalCuts &external_cuts() const { return ext_cuts; }

    const LP::CutMonitor &cut_monitor() const { return cut_mon; }

    const LP::ActiveTour &get_active_tour() const { return active_tour; }

    /// Average number of iterations per primal_pivot.
    int avg_itcount() const { return sum_it_count / num_nd_pivots; }

    double active_tourlen() const { return active_tour.length(); }


    friend class CMR::Solver;

private:
    void handle_aug_pivot(std::vector<int> tour_nodes, Basis aug_base);
    void set_active_tour(std::vector<int> tour_nodes);

    void prune_slacks();

    void rebuild_basis();

    void purge_gmi();

    Graph::CoreGraph &core_graph;
    Data::SupportGroup supp_data;

    Sep::ExternalCuts ext_cuts;

    LP::CutMonitor cut_mon;

    ActiveTour active_tour;

    std::vector<double> lp_edges;
    std::vector<double> pi_vals;
    std::vector<double> feas_stat;

    int num_nd_pivots = 0;
    int sum_it_count = 0;

    int prev_numrows;
};


}
}

#endif
