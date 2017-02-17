/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Separating over a pool of cuts.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_POOL_SEP_H
#define CMR_POOL_SEP_H

#include "hypergraph.hpp"
#include "datagroups.hpp"
#include "process_cuts.hpp"

#include <unordered_map>
#include <vector>

namespace CMR {
namespace Sep {

/// Class for cut pool separation.
class PoolCuts {
public:
    PoolCuts(ExternalCuts &EC_,
             const std::vector<Graph::Edge> &core_edges_,
             Data::SupportGroup &s_dat);

    bool find_cuts(); //!< Search cut_pool for cuts violated by supp_data.
    bool price_cuts(); //!< Find cut_slacks; return true if some are negative.

    const CutQueue<HyperGraph> &pool_q() const { return hg_q; }
    const std::vector<double> &get_slacks() const { return cut_slacks; }

private:
    ExternalCuts &EC;
    const std::vector<Graph::Edge> &core_edges;
    Data::SupportGroup &supp_data;

    /// If cuts are found they are removed from the cut pool and placed here.
    CutQueue<HyperGraph> hg_q;

    /// The slack on each cut relative to the solution in supp_data.
    std::vector<double> cut_slacks;

    /// Cut values associated to cliques.
    std::unordered_map<Clique, double> clique_vals;

    void price_cliques(); //!< Populate clique_vals.
};
}
}

#endif
