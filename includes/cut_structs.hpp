#ifndef CMR_CUT_STRUCTS_H
#define CMR_CUT_STRUCTS_H

#include "tooth.hpp"

#include <vector>

namespace CMR {
namespace Sep {

/** Structure for storing blossom inequalities from exact primal separation.
 * This structure stores blossom inequalities found by the exact primal
 * separation algorithm of Letchford-Lodi in Primal Separation Algorithms
 * Revisited. 
 */
struct blossom {
    blossom(std::vector<int> &_handle, int _cut_edge, double _val) :
        handle(_handle), cut_edge(_cut_edge), cut_val(_val){}

  
    std::vector<int> handle;
  
    int cut_edge;

    double cut_val;
};

/** Structure for storing simple DP inequalities.
 * In all cases, all numbers and indices refer to position in some `tour_nodes`
 * vector. They must be translated by deferencing, e.g., if `(i, j)` is an
 * edge in \p nonneg_edges, it refers to the edge 
 * `tour_nodes[i], tour_nodes[j]`.
 */
struct dominoparity {
    dominoparity() = default;
    dominoparity(std::vector<CMR::SimpleTooth> &_used_teeth,
                 std::vector<int> &_degree_nodes,
                 std::vector<IntPair> &_nonneg_edges) :
        used_teeth(_used_teeth), degree_nodes(_degree_nodes),
        nonneg_edges(_nonneg_edges) {}

    /** Simple tooth inequalities to be aggregated to get the simple DP ineq. */
    std::vector<CMR::SimpleTooth> used_teeth;

    /** The handle of the inequality, nodes where the degree eqns are used. */
    std::vector<int> degree_nodes;

    /** Edges \f$ e \f$ for which \f$ x_e \ge 0 \f$ is used. */
    std::vector<IntPair> nonneg_edges;
};

}  
}



#endif
