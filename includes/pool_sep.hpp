/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Separating over a pool of cuts.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_POOL_SEP_H
#define CMR_POOL_SEP_H

#include "active_tour.hpp"
#include "hypergraph.hpp"
#include "datagroups.hpp"
#include "process_cuts.hpp"
#include "cc_lpcuts.hpp"

#include <unordered_map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

/// Class for cut pool separation.
class PoolCuts {
public:
    PoolCuts(ExternalCuts &EC_,
             const std::vector<Graph::Edge> &core_edges_,
             const LP::ActiveTour &active_tour_,
             Data::SupportGroup &s_dat);

    /// Search pool for primal violated cuts, enqueueing them if found.
    bool find_cuts();

    bool tighten_pool(); //!< Attempt to tighten cuts in the pool.

    /// Populate lp_slacks and tour_slacks; return true if primal cuts found.
    bool price_cuts(bool tighten);

    CutQueue<HyperGraph> &pool_q() { return hg_q; }

    LPcutList &tighten_q() { return tight_q; }

    const std::vector<double> &get_lp_slacks() const { return lp_slacks; }
    const std::vector<double> &get_tour_slacks() const { return tour_slacks; }

    bool filter_primal = true; //!< Should only tight cuts be included.
    int verbose = 0;

private:
    void price_cliques(); //!< Populate clique_vals.

    /// Are we interested in a cut based on its tour/lp slack.
    bool slack_of_interest(double lp_slack, double tour_slack,
                           bool tighten);

    ExternalCuts &EC;
    const std::vector<Graph::Edge> &core_edges;
    const LP::ActiveTour &active_tour;
    const std::vector<double> &tour_edges;
    Data::SupportGroup &supp_data;
    Graph::AdjList tour_adj;

    /// If cuts are found they are removed from the cut pool and placed here.
    CutQueue<HyperGraph> hg_q;

    /// Successfully tightened cuts are placed here.
    LPcutList tight_q;

    /// The slack on each cut relative to the solution in supp_data.
    std::vector<double> lp_slacks;

    /// The slack on each cut relative to the tour in tour_edges.
    std::vector<double> tour_slacks;

    /// Cut values associated to cliques for lp (first) and tour (second).
    std::unordered_map<Clique, std::pair<double, double>> clique_vals;


};
}
}

#endif
