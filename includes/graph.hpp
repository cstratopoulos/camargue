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

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace CMR {

/// Classes and functions for working with graphs.
namespace Graph {

struct Edge : EndPts {
    Edge() : removable(false) {}
    Edge(int e0, int e1, int _len) :
        EndPts(e0, e1), len(_len), removable(false) {}

    bool operator==(const Edge &rhs) const {
        return ((end[0] == rhs.end[0]) && (end[1] == rhs.end[1]));
    }

    int len;
    bool removable;
};

/** Wrapper to the Concorde CCtsp_lpgraph structure.
 * This class constructs a CCtsp_lpgraph which corresponds to a tour specified
 * by the constructor arguments. It is used to check whether cuts found by 
 * Concorde standard heuristics are tight at the current tour. 
 */
class TourGraph {
public:
    TourGraph() noexcept;
    
    TourGraph(const std::vector<int> &tour_edges,
              const std::vector<Edge> &edges,
              const std::vector<int> &perm);

    TourGraph(TourGraph &&T) noexcept;
    TourGraph &operator=(TourGraph &&T) noexcept;
    
    ~TourGraph();

    TourGraph(const TourGraph &T) = delete;
    TourGraph &operator=(const TourGraph &T) = delete;

    

    CCtsp_lpgraph* pass_ptr() { return &L; }
    double* tour_array() { return &d_tour[0]; }
    int node_count() const { return L.ncount; }
  
private:
  std::vector<double> d_tour;
  CCtsp_lpgraph L;
};

struct AdjObj {
    AdjObj() = default;
    AdjObj(int otherend, int index, double _val) :
        other_end(otherend), edge_index(index), val(_val) {}
    
    int other_end;
    int edge_index;
    double val;

    bool operator==(const AdjObj &rhs) const
    {
        return (other_end == rhs.other_end
                && edge_index == rhs.edge_index
                && val == rhs.val);
    }
};

struct Node {
    Node() : mark(0) {}
    
    int degree() const { return neighbors.size(); }
    std::vector<AdjObj> neighbors;
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

    bool connected(std::vector<int> &island, int start_node);

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

/// Graph object representing the edges in a core lp relaxation.
class CoreGraph {
public:
    CoreGraph() = default;

    /// Construct a CoreGraph from a length func and a c-array node-node elist.
    CoreGraph(int ncount, int ecount, const int *elist,
              const std::function<double(int, int)> edgelen);

    /// Construct a CoreGraph containing the nodes of a TSP tour.
    CoreGraph(const std::vector<int> &tour_nodes,
              const std::function<double(int, int)> edgelen);

    int node_count() const { return nodecount; }
    int edge_count() const { return edges.size(); }

    int find_edge_ind(int end0, int end1) const;

    Edge get_edge(int index) const { return edges[index]; }
    
    std::vector<Edge> &get_edges() { return edges; }
    const std::vector<Edge> &get_edges() const { return edges; }

    void get_elist(std::vector<int> &elist, std::vector<int> &elen) const;

    const AdjList &get_adj() const { return adj_list; }
    
    void add_edge(int end0, int end1, int len);
    void add_edge(const Edge &e);

private:
    std::vector<Edge> edges;
    AdjList adj_list;
    int nodecount;
};

//TODO get rid of all the versions that need island/deltacount, etc. 
std::vector<int> delta_inds(const std::vector<int> &node_list,
                            const std::vector<Edge> &edges,
                            int ncount);

void get_delta (const std::vector<int> &nodelist, std::vector<Edge> &elist,
		int *deltacount_p, std::vector<int> &delta,
		std::vector<int> &marks);

void get_delta (int nsize, int *nlist, int ecount, int *elist,
		int *deltacount, int *delta, int *node_marks);

void get_delta (const int interval_start, const int interval_end,
		const std::vector<int> &tour_nodes,
		const std::vector<int> &elist,
		int &deltacount,  std::vector<int> &delta,
		std::vector<int> &node_marks);



///////////////// TEMPLATE METHOD IMPLEMENTATIONS /////////////////////////////

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
