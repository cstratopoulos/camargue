#ifndef CMR_PRICER_H
#define CMR_PRICER_H

#include "core_lp.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "price_util.hpp"
#include "fixed64.hpp"
#include "edgehash.hpp"

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

#include <queue>
#include <unordered_map>
#include <vector>

namespace CMR {

/// Matters related to pricing sets of edges.
namespace Price {

/// Get reduced costs for edges not in the core lp. 
class Pricer {
public:
    
    Pricer(LP::CoreLP &core, const Data::Instance &_inst,
           Data::GraphGroup &graphgroup); //!< Construct a Pricer.

    Pricer(const Pricer &P) = delete; //!< Deleted copy constructor.
    Pricer &operator=(const Pricer &P) = delete; //!< Deleted copy assign.

    ~Pricer(); //!< Destruct and free resource handles.
    
    ScanStat gen_edges(LP::PivType piv_stat); //!< Generate/add edges to core.

    void exact_lb();

    void price_edges(std::vector<PrEdge> &target_edges,
                     bool compute_duals); //!< Price a list of edges.

    int queue_size() const { return edge_q.size(); } //<! Size of edge queue.

private:
    void sort_q(); //!< Sort the queue of edges by minimum reduced costs.
    std::vector<Graph::Edge> get_pool_chunk(); //!< Get at most AddBatch edges.
    
    LP::CoreLP &core_lp; //!< The LP relaxation to query/modify.
    const Data::Instance &inst; //!< To get lengths for edges not in core_lp.
    const Sep::ExternalCuts &ext_cuts;  //!< For computing duals.

    Data::GraphGroup &graph_group; //!< Graph data for the core_lp.

    const int gen_max; //!< The max number of edges to generate at a time.
    
    std::vector<int> gen_elist; //<! Raw node-node list of generated edges.
    std::vector<int> gen_elen;  //<! Unused dummy parameter to pass.

    std::vector<double> node_pi; //!< pi values for degree eqns.
    std::vector<double> node_pi_est; //!< estimated node pi for dominos.
    
    std::vector<double> cut_pi; //!< pi values for cuts.

    /// Dual values/multiplicities for Cliques. 
    std::unordered_map<Sep::Clique, double> clique_pi;

    CCtsp_edgegenerator eg_inside; //<! Concorde 50-nearest edge generator.
    CCtsp_edgegenerator eg_full; //<! Concorde complete graph edge generator.

    util::EdgeHash edge_hash; //<! Hash table for tracking generated edges.

    std::vector<PrEdge> price_elist; //!< Compute exact RCs for these edges.
    std::vector<PrEdge> edge_q; //!< Queue of edges for inclusion.

};

}
}

#endif
