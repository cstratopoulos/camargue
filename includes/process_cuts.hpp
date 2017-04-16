/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Structures for storing and processing cuts.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_PROCESS_CUTS_H
#define CMR_PROCESS_CUTS_H

#include "datagroups.hpp"
#include "cc_lpcuts.hpp"
#include "graph.hpp"
#include "cut_structs.hpp"
#include "lp_util.hpp"
#include "util.hpp"

#include <memory>
#include <vector>
#include <limits>
#include <utility>

namespace CMR {
namespace Sep {

/// SparseRow corresponding to Concorde cut.
LP::SparseRow get_row(const CCtsp_lpcut_in &cc_cut,
                      const std::vector<int> &perm,
                      const Graph::CoreGraph &core_graph);

/// SparseRow corresponding to simple DP inequality.
LP::SparseRow get_row(const dominoparity &dp_cut,
                      const std::vector<int> &tour_nodes,
                      const Graph::CoreGraph &core_graph);

/// SparseRow corresponding to the handle and tooth edges of a blossom.
LP::SparseRow get_row(const std::vector<int> &handle_delta,
                      const std::vector<std::vector<int>> &tooth_edges,
                      const Graph::CoreGraph &core_graph);

/// Function template for determining the activity or lhs of a vector on a row.
template<typename number_type>
double get_activity(const std::vector<number_type> &x,
                    const LP::SparseRow &R)
{
    double result = 0.0;
    for(auto i = 0; i < R.rmatind.size(); i++){
        int index = R.rmatind[i];
        result += x[index] * R.rmatval[i];
    }
    return result;
}


/// Gets the indices of the teeth for an ex_blossom \p B relative to \p edges.
std::vector<int> teeth_inds(const ex_blossom &B,
                            const std::vector<double> &tour_edges,
                            const std::vector<double> &lp_vec,
                            const std::vector<Graph::Edge> &edges,
                            int ncount);

/// Like the other version, but if we already have handle_delta.
std::vector<int> teeth_inds(const ex_blossom &B,
                            const std::vector<double> &tour_edges,
                            const std::vector<double> &lp_vec,
                            const std::vector<Graph::Edge> &edges,
                            int ncount, const std::vector<int> &handle_delta);


/// Returns true if the blossom is invalid for some reason.
bool bad_blossom(const ex_blossom &B,
                 const std::vector<double> &tour_edges,
                 const std::vector<double> &lp_vec,
                 const std::vector<Graph::Edge> &edges, int ncount);

}
}



#endif
