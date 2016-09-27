#ifndef PSEP_GRAPH_H
#define PSEP_GRAPH_H

#include <vector>
#include<queue>
#include <iostream>
#include<memory>

#include<math.h>

#include "PSEP_util.hpp"

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

namespace PSEP {
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

    operator bool() const { return node_count > 0; }
    
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

    //TODO: This should just be a member function of SupportGraph
    int build_s_graph (int node_count, int edge_count,
		       std::vector<Edge> &edges,
		       std::vector<int> &support_indices,
		       std::vector<double> &m_lp_edges,
		       SupportGraph *G_s);
  }

  struct CC {
    struct GH {    
      static void grab_cut_dfs(CC_GHnode *n, std::vector<int> &cut_nlist){
	for(int i = 0; i < n->listcount; i++)
	  cut_nlist.push_back(n->nlist[i]);

	for(n = n->child; n; n = n->sibling)
	  grab_cut_dfs(n, cut_nlist);
      }
    
      static void odd_cut_dfs(CC_GHnode *n, double *min_val_p,
			      CC_GHnode **best_n){
	if(n->parent){
	  if(n->cutval < *min_val_p && (n->ndescendants % 2) == 1){
	    *best_n = n;
	    *min_val_p = n->cutval;
	  }
	}

	if(*min_val_p >= PSEP::LP::EPSILON)
	  for(n = n->child; n; n = n->sibling)
	    odd_cut_dfs(n, min_val_p, best_n);
      }

      static void get_odd_cut(CC_GHtree *T, std::vector<int> &cut_nlist){
	double min_val = 1 - PSEP::LP::EPSILON;
	CC_GHnode *best_node = (CC_GHnode *) NULL;

	if(T && T->root)
	  odd_cut_dfs(T->root, &min_val, &best_node);

	if((min_val < 1 - PSEP::LP::EPSILON) && best_node)
	  grab_cut_dfs(best_node, cut_nlist);
      }
    };
  };
}

#endif
