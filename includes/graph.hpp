/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Header for classes/structures/functions to work with graphs.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_GRAPH_H
#define CMR_GRAPH_H

#include "util.hpp"

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cmath>

namespace CMR {

/// Classes and functions for working with graphs.
namespace Graph {

/// Representing graph edges and costs.
struct Edge : EndPts {
    Edge() : removable(false) {}

    /// Construct by endpoints and length.
    Edge(int e0, int e1, int _len) :
        EndPts(e0, e1), len(_len), removable(false) {}

    bool operator==(const Edge &rhs) const {
        return ((end[0] == rhs.end[0]) && (end[1] == rhs.end[1]));
    }

    int len; //!< The edge length.
    bool removable; //!< For use by certain erase-remove procedures.
};

/// Get a node-node elist representation of a list of edges.
void get_elist(const std::vector<Edge> &edges, std::vector<int> &elist,
               std::vector<int> &ecap);

/// Object used to represent adjacency in a Graph::AdjList.
/// For a fixed node `x`, an AdjObj is used to represent adjacency between
/// `x` and another node `y` in the graph.
struct AdjObj {
    AdjObj() = default;

    AdjObj(int otherend, int index, double _val) :
        other_end(otherend), edge_index(index), val(_val) {}

    int other_end; //!< The node `y` in the notation above.
    int edge_index; //!< The index wrt a fixed Graph::Edge vector of `(x, y)`.
    double val; //!< Can represent length or some other value on the edge.

    bool operator==(const AdjObj &rhs) const
    {
        return (other_end == rhs.other_end
                && edge_index == rhs.edge_index
                && val == rhs.val);
    }
};

/// A vertex in a Graph::AdjList graph.
/// For a node \f$ x \f$, this effectively stores \f$ \delta(x) \f$.
struct Node {
    Node() : mark(0) {}

    int degree() const { return neighbors.size(); }
    std::vector<AdjObj> neighbors; //!<
    int mark;
};

/// Representation of a graph as an adjacency list.
struct AdjList {
    AdjList() = default;

    /// An AdjList with \p ncount nodes for all the edges in \p ref_elist.
    AdjList(int ncount, const std::vector<Edge> &ref_elist);

    /// An AdjList for a vector of structs derived from EndPt.
    template <typename EndPt_type>
    AdjList(int ncount, const std::vector<EndPt_type> &elist);

    /// A support graph type AdjList.
    AdjList(int ncount,
            const std::vector<Edge> &ref_elist,
            const std::vector<double> &ref_elist_caps,
            const std::vector<int> &keep_indices);

    AdjList(AdjList &&AL) noexcept;

    AdjList &operator=(AdjList &&AL) noexcept;

    /// Is the graph connected.
    bool connected(std::vector<int> &island, int start_node);

    /// Performs a depth-first beginning with \p start_node.
    void dfs(int start_node, std::vector<int> &island);

    /// Get a pointer to the AdjObj with end points \p end0 and \p end1.
    /// @returns `nullptr` if not found, else a pointer to the AdjObj.
    const AdjObj *find_edge(int end0, int end1) const
    {
        for (const AdjObj &a : nodelist[end0].neighbors)
            if (a.other_end == end1)
                return &a;
        return nullptr;
    }

    /// Add the edge with end points \p end0 \p end1 to the AdjList.
    void add_edge(int end0, int end1, int index, double val);

    int node_count;
    int edge_count;

    std::vector<Node> nodelist;
};


/** @name Functions for getting cuts in a graph.
 * These functions take as input a graph \f$ G \f$ and a node set \f$ X \f$,
 * returning a representation of \f$ \delta(X) \f$ in \f$ G \f$. The returned
 * vector \p delta_inds will store indices of edges in the cut wrt the input
 * graph.
 */
///@{

/// Cut set \p node_list, with graph specified by \p edges with \p ncount nodes.
std::vector<int> delta_inds(const std::vector<int> &node_list,
                            const std::vector<Edge> &edges,
                            int ncount);

/// As above but with node-node list \p elist representing graph edges.
std::vector<int> delta_inds(const std::vector<int> &node_list,
                            const std::vector<int> &elist,
                            int ncount);

///@}



///////////////// TEMPLATE METHOD IMPLEMENTATIONS /////////////////////////////

/**
 * @tparam EndPt_type the edge representation being used. Should be derived
 * from CMR::EndPt.
 */
template <typename EndPt_type>
AdjList::AdjList(int ncount, const std::vector<EndPt_type> &elist) try
    : node_count(ncount), edge_count(elist.size()),
      nodelist(std::vector<Node>(node_count))
{
    for (int i = 0; i < edge_count; ++i) {
        const auto &e = elist[i];

        nodelist[e.end[0]].neighbors.emplace_back(e.end[1], i, 0.0);
        nodelist[e.end[1]].neighbors.emplace_back(e.end[0], i, 0.0);
    }
} catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    throw std::runtime_error("AdjList EndPt_type constructor failed.");
}

}
}

#endif
