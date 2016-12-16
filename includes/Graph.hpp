#ifndef CMR_GRAPH_H
#define CMR_GRAPH_H

#include <vector>
#include <queue>
#include <iostream>
#include <memory>

#include <math.h>

#include "util.hpp"

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
  
  int end[2];
  int len;
  bool removable;
};

struct Graph {
  Graph() : node_count(0), edge_count(0) {}
    
  void print_edges();

  int node_count;
  int edge_count;
  std::vector<Edge> edges;
  IntPairMap edge_lookup;

  void print_edge(int i) const {
    if(i >= edge_count)
      std::cout << "Edge out of range" << std::endl;
    else
      std::cout << "Edge " << i << ": [" << edges[i].end[0]
		<< ", " << edges[i].end[1] << "]" << std::endl;
  };
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
}

#endif
