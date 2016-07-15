#ifndef PSEP_GRAPH_H
#define PSEP_GRAPH_H

#include <vector>
#include<queue>
#include<utility>
#include <iostream>

extern "C" {
#include "../programs/concorde/concorde.h"
}

#include "PSEP_util.h"


struct Edge {
  Edge() {}
  Edge(int e0, int e1, int _len);
  
  int end[2];
  int len;

  bool operator<(const Edge& val) const {
    return len < val.len;
  }

  static bool ptr_compare(Edge *e0, Edge *e1);
};

struct Graph {
Graph() : node_count(0), edge_count(0) {}

  void print_edges();

  int node_count;
  int edge_count;
  std::vector<Edge> edges;

  void print_edge(int i){
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


struct G_Utils {
  static int connected (SupportGraph *G, int *icount,
			std::vector<int> &island, int starting_node);
  static void dfs (int n, SupportGraph *G, int *icount,
		   std::vector<int> &island);  
  static void get_delta (std::vector<int> &nodelist, std::vector<Edge> &elist,
		  int *deltacount_p, std::vector<int> &delta,
		  std::vector<int> &marks);
  static void get_delta (int nsize, int *nlist, int ecount, int *elist,
			 int *deltacount, int *delta, int *edge_marks);

  static int build_s_graph (int node_count, int edge_count,
			    std::vector<Edge> &edges,
			    std::vector<int> &support_indices,
			    std::vector<double> &m_lp_edges,
			    SupportGraph *G_s);
};

struct CC {
  struct GH {
    class Comparator {
    public:
      bool operator() (const CC_GHnode *a, const CC_GHnode *b){
	return a->cutval < b->cutval;
      }
    };

    typedef std::priority_queue<CC_GHnode*, std::vector<CC_GHnode*>,
      Comparator> cut_pq;

    static void grab_cut_dfs(CC_GHnode *n, std::vector<int> &cut_nlist);
    static void pq_dfs(CC_GHnode *n, const int max_cutcount, cut_pq &pq);
    static void get_all_toothlists(CC_GHtree *T, const int max_cutcount,
				   cut_pq &pq,
				   std::vector<std::vector<int> >
				   &toothlists);
  };
};

#endif
