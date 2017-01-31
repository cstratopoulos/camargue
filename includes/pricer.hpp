#ifndef CMR_PRICER_H
#define CMR_PRICER_H

#include "core_lp.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "util.hpp"
#include "price_util.hpp"
#include "err_util.hpp"
#include "fixed64.hpp"
#include "edgehash.hpp"

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

#include <memory>
#include <queue>
#include <stdexcept>
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

    Pricer(const Pricer &P) = delete;
    Pricer &operator=(const Pricer &P) = delete;

    ~Pricer(); //!< Destruct and free resource handles.
    
    ScanStat gen_edges(LP::PivType piv_stat); //!< Generate/add edges to core.

    //void exact_lb();
    template <typename numtype>
    void price_edges(std::vector<PrEdge<numtype>> &target_edges,
                     std::unique_ptr<LP::DualGroup<numtype>> &duals);

private:
    std::vector<Graph::Edge> pool_chunk(std::vector<PrEdge<double>> &edge_q);
    
    // bool f64_gen_edges(const std::vector<util::Fixed64> &node_pi_est,
    //                    std::vector<PrEdge<double>64> &gen_edges,
    //                    int &loop1, int &loop2);

    // void f64_price_edges(std::vector<PrEdge<double>64> &target_edges,
    //                      std::vector<util::Fixed64> &node_pi,
    //                      std::vector<util::Fixed64> &node_pi_est,
    //                      std::vector<util::Fixed64> &cut_pi,
    //                      std::unordered_map<Sep::Clique,
    //                      util::Fixed64> &clique_pi);
    
    LP::CoreLP &core_lp; //!< The LP relaxation to query/modify.
    const Data::Instance &inst; //!< To get lengths for edges not in core_lp.
    const Sep::ExternalCuts &ext_cuts;  //!< For computing duals.

    Data::GraphGroup &graph_group; //!< Graph data for the core_lp.

    const int gen_max; //!< The max number of edges to generate at a time.
    
    std::vector<int> gen_elist; //<! Raw node-node list of generated edges.
    std::vector<int> gen_elen;  //<! Unused dummy parameter to pass.

    std::unique_ptr<LP::DualGroup<double>> reg_duals;
    std::unique_ptr<LP::DualGroup<util::Fixed64>> ex_duals;

    CCtsp_edgegenerator eg_inside; //<! Concorde 50-nearest edge generator.
    CCtsp_edgegenerator eg_full; //<! Concorde complete graph edge generator.

    util::EdgeHash edge_hash; //<! Hash table for tracking generated edges.
};

//////////////////// TEMPLATE METHOD IMPLEMENTATIONS //////////////////////////

template <typename numtype>
void Pricer::price_edges(std::vector<PrEdge<numtype>> &target_edges,
                         std::unique_ptr<LP::DualGroup<numtype>> &duals)
{
    using std::vector;
    using std::unique_ptr;
    using Dual = LP::DualGroup<numtype>;
    using CutType = Sep::HyperGraph::Type;

    std::runtime_error err("Problem in Pricer::price_edges.");

    if (!duals)
        try {
            duals = util::make_unique<Dual>(false, core_lp,
                                            core_lp.external_cuts());
        } CMR_CATCH_PRINT_THROW("getting duals", err);

    vector<numtype> &node_pi = duals->node_pi;
    vector<numtype> &cut_pi = duals->cut_pi;
    std::unordered_map<Sep::Clique, numtype> &clique_pi = duals->clique_pi;

    for (auto &e : target_edges)
        e.redcost = inst.edgelen(e.end[0], e.end[1]) - node_pi[e.end[0]]
        - node_pi[e.end[1]];

    Graph::AdjList price_adjlist;

    try  {
        price_adjlist = Graph::AdjList(inst.node_count(), target_edges);
    } CMR_CATCH_PRINT_THROW("Couldn't build price adjlist.", err);

    vector<Graph::Node> &price_nodelist = price_adjlist.nodelist;
    
    const std::vector<int> &def_tour = ext_cuts.get_cbank().ref_tour();
    int marker = 0;
    
    for (const std::pair<Sep::Clique, numtype> &kv : clique_pi) {
        const Sep::Clique &clq = kv.first;
        numtype pival = kv.second;

        if (pival != 0.0) {
            numtype add_back = pival + pival;
            ++marker;
            
            for (int j : clq.node_list(def_tour)) {
                for (Graph::AdjObj &nbr : price_nodelist[j].neighbors)
                    if (price_nodelist[nbr.other_end].mark == marker)
                        target_edges[nbr.edge_index].redcost += add_back;
                    
                price_nodelist[j].mark = marker;
            }
        }
    }

    const vector<Sep::HyperGraph> &cutlist = ext_cuts.get_cuts();
    vector<int> rmatind;
    vector<double> rmatval;

    for (int i = 0; i < cutlist.size(); ++i) {
        numtype pival = cut_pi[i];
        const Sep::HyperGraph &H = cutlist[i];
        
        if (H.cut_type() == CutType::Non)
            throw std::logic_error("Called pricing w Non HyperGraph present.");
        
        if (H.cut_type() != CutType::Domino)
            continue;

        if (pival == 0)
            continue;
        
        try {
            H.get_coeffs(target_edges, rmatind, rmatval);
        } CMR_CATCH_PRINT_THROW("geting domino price edge coeffs", err);

        for (int j = 0; j < rmatind.size(); ++j)
            util::add_mult(target_edges[rmatind[j]].redcost,
                           pival, -rmatval[j]);
    }
}

}
}

#endif
