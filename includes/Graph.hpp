#ifndef CMR_GRAPH_H
#define CMR_GRAPH_H

#include "price_util.hpp"
#include "util.hpp"

#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace CMR {


struct Edge {
    Edge() : removable(false) {}
    Edge(int e0, int e1, int _len);

    bool operator==(const Edge &rhs) const {
        return ((end[0] == rhs.end[0]) && (end[1] == rhs.end[1]));
    }

    std::array<int, 2> end;
    int len;
    bool removable;
};

struct s_adjobj {
  int other_end;
  int edge_index;
  double lp_weight;
};

struct SNode {
  int s_degree;
  s_adjobj *adj_objs;
  int mark;
};

struct SupportGraph {
  SupportGraph() :
    node_count(0), edge_count(0), nodelist((SNode *) NULL),
    adjspace((s_adjobj *) NULL) {}

  ~SupportGraph(){
    if(nodelist) free(nodelist);
    if(adjspace) free(adjspace); 
  }

  
  int node_count;
  int edge_count;
  SNode *nodelist;
  s_adjobj *adjspace;
};

/** Wrapper to the Concorde CCtsp_lpgraph structure.
 * This class constructs a CCtsp_lpgraph which corresponds to a tour specified
 * by the constructor arguments. It is used to check whether cuts found by 
 * Concorde standard heuristics are tight at the current tour. 
 */
class TourGraph {
public:
    TourGraph(const std::vector<int> &tour_edges,
              const std::vector<CMR::Edge> &edges,
              const std::vector<int> &perm);
    ~TourGraph();

    CCtsp_lpgraph* pass_ptr() { return &L; }
    double* tour_array() { return &d_tour[0]; }
    int node_count() const { return L.ncount; }
  
private:
  CCtsp_lpgraph L;
  std::vector<double> d_tour;
};


namespace GraphUtils {
int connected (SupportGraph *G, int *icount,
	       std::vector<int> &island, int starting_node);
void dfs (int n, SupportGraph *G, int *icount,
	  std::vector<int> &island);  
void get_delta (const std::vector<int> &nodelist, std::vector<Edge> &elist,
		int *deltacount_p, std::vector<int> &delta,
		std::vector<int> &marks);
void get_delta (int nsize, int *nlist, int ecount, int *elist,
		int *deltacount, int *delta, int *edge_marks);

void get_delta (const int interval_start, const int interval_end,
		const std::vector<int> &tour_nodes,
		const std::vector<int> &elist,
		int &deltacount,  std::vector<int> &delta,
		std::vector<int> &edge_marks);

//TODO: This should just be a member function of SupportGraph
int build_s_graph (int node_count,
		   const std::vector<Edge> &edges,
		   const std::vector<int> &support_indices,
		   const std::vector<double> &m_lp_edges,
		   SupportGraph *G_s);
}

namespace GraphUtils {

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
    Node() = default;
    Node(int degree_estimate) { neighbors.reserve(degree_estimate); }
    
    int degree() const { return neighbors.size(); }
    std::vector<AdjObj> neighbors;
    int mark;
};

struct AdjList {
    AdjList() = default;
    
    AdjList(int ncount, const std::vector<CMR::Edge> &ref_elist);

    AdjList(int ncount, const std::vector<Price::edge> &price_elist);
    
    AdjList(int ncount,
            const std::vector<CMR::Edge> &ref_elist,
            const std::vector<double> &ref_elist_caps,
            const std::vector<int> &keep_indices);

    const AdjObj *find_edge(int end0, int end1) const
    {
        for (const AdjObj &a : nodelist[end0].neighbors)
            if (a.other_end == end1)
                return &a;
        return nullptr;
    }

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

    CMR::Edge &get_edge(int index) { return edges[index]; }
    
    std::vector<CMR::Edge> &get_edges() { return edges; }
    const std::vector<CMR::Edge> &get_edges() const { return edges; }
    
    void add_edge(int end0, int end1, int len);    

private:
    std::vector<CMR::Edge> edges;
    AdjList adj_list;
    int nodecount;
};

}
}

#endif
