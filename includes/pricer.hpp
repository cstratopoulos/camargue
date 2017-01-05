#ifndef CMR_PRICER_H
#define CMR_PRICER_H

#include "lp_interface.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "price_util.hpp"

#include <queue>
#include <unordered_map>
#include <vector>

namespace CMR {

/** Namespace for manners related to pricing sets of edges. */
namespace Price {

constexpr int nearest_factor = 50;

/** Class for pricing edges not in the CoreLP. */
class Pricer {
public:
    /** Construct a Pricer.
     * @param[in] _relax the Relaxation for grabbing dual values
     * @param[in] _inst the TSP instance for generating edges
     * @param[in] _ext_cuts the HyperGraph representation of the cuts in 
     * \p _relax.
     */
    Pricer(const LP::Relaxation &_relax, const Data::Instance &_inst,
           const Sep::ExternalCuts &_ext_cuts);

    /** Generate edges to add to the CoreLP.
     * @param[in] piv_stat the PivType from the solution when edge generation 
     * is called. Used to determine which pricing protocol will be used.
     * @returns a ScanStat indicating the outcome of the pricing scan.
     */
    ScanStat add_edges(LP::PivType piv_stat);

private:
    /// Scan the 50-nearest edges.
    /// @returns true iff negative reduced cost edges were found.
    ScanStat partial_scan();

    /// Scan the full edge set.
    /// @returns true iff negative reduced cost edges were found.
    ScanStat full_scan(); 
    
    const LP::Relaxation &relax;
    const Data::Instance &inst;
    const Sep::ExternalCuts &ext_cuts;

    std::vector<double> node_pi;
    std::vector<double> node_pi_est;
    
    std::vector<double> cut_pi;

    std::unordered_map<Sep::Clique, double> clique_pi;
    
    int last_ind;

    util::c_array_ptr nearest_elist;
    int nearest_ecount;


    std::vector<edge> edge_q;

};

}
}

#endif
