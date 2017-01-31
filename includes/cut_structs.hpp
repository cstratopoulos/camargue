#ifndef CMR_CUT_STRUCTS_H
#define CMR_CUT_STRUCTS_H

#include "util.hpp"
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

struct ToothBody : Segment {
    ToothBody() = default;
    ToothBody(int _start, int _end, double _slack) :
        Segment(_start, _end), slack(_slack) {}
    ToothBody(int _start, int _end) : Segment(_start, _end), slack(1.0) {}

    double slack;
};

struct SimpleTooth {
    SimpleTooth(int _root, int _body_start, int _body_end, double _slack) :
        root(_root), body_start(_body_start), body_end(_body_end),
        slack(_slack) {}

    SimpleTooth(int _root, ToothBody &seg, double _slack) :
        root(_root), body_start(seg.start), body_end(seg.end), slack(_slack) {}

    using Ptr = std::unique_ptr<SimpleTooth>;

    int root, body_start, body_end;
    int cutgraph_index;
    double slack;

    enum Type {
        LeftAdj = 0,
        RightAdj = 1,
        Dist = 2
    };

    Type type() const
        {
            if (body_start == root + 1)
                return LeftAdj;
            if (body_end + 1 == root)
                return RightAdj;
            return Dist;
        }

    int body_size() const { return body_end - body_start + 1; }
    bool body_contains(int i) const { return body_start <= i && i <= body_end; }
    bool is_subset_of(const SimpleTooth &T) const {
        return root == T.root &&
        T.body_start <= body_start &&
        body_end <= T.body_end;
    }
    
};

/** Structure for storing simple DP inequalities.
 * In all cases, all numbers and indices refer to position in some `tour_nodes`
 * vector. They must be translated by deferencing, e.g., if `(i, j)` is an
 * edge in \p nonneg_edges, it refers to the edge 
 * `tour_nodes[i], tour_nodes[j]`.
 */
struct dominoparity {
    dominoparity() = default;
    dominoparity(std::vector<SimpleTooth> &_used_teeth,
                 std::vector<int> &_degree_nodes,
                 std::vector<IntPair> &_nonneg_edges) :
        used_teeth(_used_teeth), degree_nodes(_degree_nodes),
        nonneg_edges(_nonneg_edges) {}

    /** Simple tooth inequalities to be aggregated to get the simple DP ineq. */
    std::vector<SimpleTooth> used_teeth;

    /** The handle of the inequality, nodes where the degree eqns are used. */
    std::vector<int> degree_nodes;

    /** Edges \f$ e \f$ for which \f$ x_e \ge 0 \f$ is used. */
    std::vector<IntPair> nonneg_edges;
};

/// Simple struct representing sparse matrix row for passing to LP solver. 
struct SparseRow {    
    std::vector<int> rmatind;
    std::vector<double> rmatval;
    char sense;
    double rhs;
    double lp_viol;
};

}  
}



#endif
