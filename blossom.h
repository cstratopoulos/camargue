#ifndef PSEP_BLOSSOM_H
#define PSEP_BLOSSOM_H

#include <vector>
#include <queue>
#include <math.h>

#include "Graph.h"
#include "lp.h"

extern "C" {
#include "../programs/concorde/concorde.h"
}

class PSEP_2match {
 public:
 PSEP_2match(std::vector<int> & _tour_edges, std::vector<double> & _lp_edges,
	      std::vector<int> & _sup_inds, std::vector<int> & _sup_elist,
	      std::vector<double> & _ecap, PSEPlp & _lp) :
  best_tour_edges(_tour_edges), m_lp_edges(_lp_edges),
    support_indices(_sup_inds), support_elist(_sup_elist),
    support_ecap(_ecap), m_lp(_lp) {}

  int separate(const int max_cutcount);
  bool q_empty() const {return pq.empty();};
  int q_size() const {return pq.size();};
  void pop(std::vector<int> &hnodes, int *cutedge_p, double *cutval_p){
    hnodes = pq.top().handle;
    *cutedge_p = pq.top().cut_edge;
    *cutval_p = pq.top().cut_val;

    pq.pop();
  };
  int add_cut(const int deltacount, std::vector<int> &delta, const int cutedge);


  
 private:
  std::vector<int> &best_tour_edges;
  std::vector<double> &m_lp_edges;
  std::vector<int> &support_indices;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;
  std::vector<double> cut_ecap;

  struct blossom {
  blossom(std::vector<int> _handle, int _edge, double _val) :
    handle(_handle), cut_edge(_edge), cut_val(_val){}
    blossom(){}

    bool operator< (const blossom &val) const {
      return cut_val < val.cut_val;
    }
    
    std::vector<int> handle;
    int cut_edge;
    double cut_val;
  };

  PSEPlp &m_lp;

  std::priority_queue<blossom> pq;
};

#endif
