#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include "blossom.h"
#include "segcuts.h"
#include "simpleDP.h"

class PSEP_Cutcall {
 public:
  PSEP_Cutcall(//its own objects
	       std::vector<Edge> & _edges, std::vector<int> & _delta,
	       std::vector<int> & _edge_marks,
	       //segments inits
	       SupportGraph &supgraph, std::vector<int> &nodes, PSEPlp &lp,
	       //2match inits
	       std::vector<int> & _tour_edges, std::vector<double> & _lp_edges,
	       std::vector<int> & _sup_inds, std::vector<int> & _sup_elist,
	       std::vector<double> & _ecap) :
  edges(_edges), delta(_delta), edge_marks(_edge_marks), best_tour_nodes(nodes),
    support_elist(_sup_elist), support_ecap(_ecap),
    segments(supgraph, _edge_marks, nodes, lp),
    blossoms(_tour_edges, _lp_edges, _sup_inds, _sup_elist, _ecap, lp),
    dominos(nodes, supgraph, _edge_marks, _sup_elist, _ecap) {}
	       
  
  int segment(int *num_added_p);
  int blossom(const int max_cutcount, int *num_added_p);
  int simpleDP(const int max_cutcount, int *num_added_p);

  int in_subtour_poly(bool *result_p);

 private:
  std::vector<Edge> &edges;
  std::vector<int> &delta;
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;

  std::vector<int> &support_elist;
  std::vector<int> &support_ecap;
  
  PSEP_Segcuts segments;
  PSEP_2match blossoms;
  PSEP_SimpleDP dominos;
  
};


#endif
