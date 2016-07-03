#ifndef PSEP_SEGCUTS_H
#define PSEP_SEGCUTS_H

#include <vector>
#include <queue>
#include <math.h>

#include "Graph.h"
#include "lp.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *
 *                     Segment Cut Separation Routines                        
 *
 * Functions exported to TSP_Solver via PSEP_Segcuts object:
 * 
 *   void separate()
 *      searches for violated subtour inequalities arising from segments of the
 *      current best tour, using support graph G_s. Cuts are stored as a list of
 *      pairs, see below
 *
 *  int add_cut(int deltacout, vector<int> delta)
 *      adds subtour cut corresponding to delta to lp
 *
 *
 *
 *
 * Non-reference member objects:
 *
 *    queue<pair<int, int>> segment_list
 *      -each violated segment cut found is stored as a pair of integers,
 *       indicating the start and endpoints of the tour interval
 *      -this is stored as the indices of the best tour vector, i.e., the pair
 *       (i, j) indicates that the segment is best_tour_nodes[i, ..., j]
 *      -at most 80 cuts are stored in the list at a time
 *      -dependent on existing tour_nodes/perm vector
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
 


class PSEP_Segcuts {
 public:
 PSEP_Segcuts(SupportGraph &supgraph, std::vector<int> &marks,
	      std::vector<int> &nodes, PSEPlp &lp) :
  G_s(supgraph), edge_marks(marks), best_tour_nodes(nodes), m_lp(lp){};
  
  void separate();
  bool q_empty() const {return seg_pq.empty();};
  int q_size() const {return seg_pq.size();};
  void seg_pop(int *start_p, int *end_p, double *viol_p){
    *start_p = seg_pq.top().start();
    *end_p = seg_pq.top().end();
    *viol_p = seg_pq.top().cut_viol();
    seg_pq.pop();
  }
  int add_cut(const int deltacount, std::vector<int> &delta);

 private:

  
  class seg {
  public:
  seg(int start, int end, double _slack) :
    viol(_slack){ends.first = start; ends.second = end;};
    
    int start() const {return ends.first;};
    int end() const {return ends.second;};
    double cut_viol() const {return viol;};

    bool operator <(const seg &val) const {
      return viol > val.viol;}
    
  private:
    std::pair<int, int> ends;
    double viol;
  };
  
  SupportGraph &G_s;
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;

  PSEPlp &m_lp;

  std::priority_queue<seg> seg_pq; 
};

#endif
