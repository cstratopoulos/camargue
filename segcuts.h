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
 *      searches for violated subtour inequalities arising from segments of
 *      current best tour, using support graph G_s. Cuts are stored as a list 
 *      of pairs, see below
 *
 *   int add_cut(int deltacout, vector<int> delta)
 *      adds subtour cut corresponding to delta to lp
 *
 *   void seg_pop (int *start_p, int *end_p, double *viol_p)
 *     -pops the seg on top of the priority queue; stores its parameters in
 *      the variables pointed to by function args
 *
 *
 *
 * Non-reference member objects:
 *
 *    priority_queue<seg> seg_pq is a priority queue of seg objects (see below)
 *
 *
 *   struct seg 
 *       objects used to store segment cuts
 *       -start/end are the indices of the start and endpoints of the segment
 *        a start of i indicates the first node in the segment is 
 *        best_tour_nodes[i]
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
 


class PSEP_Segcuts {
 public:
 PSEP_Segcuts(SupportGraph &supgraph, std::vector<int> &marks,
	      std::vector<int> &nodes, PSEPlp &lp) :
  G_s(supgraph), edge_marks(marks), best_tour_nodes(nodes), m_lp(lp){};
  
  void separate(const int max_cuts);
  bool q_empty() const {return pq.empty();};
  int q_size() const {return pq.size();};
  void pop(int *start_p, int *end_p, double *viol_p){
    *start_p = pq.top().start();
    *end_p = pq.top().end();
    *viol_p = pq.top().viol;
    pq.pop();
  }
  int add_cut(const int deltacount, std::vector<int> &delta);

 private:  
  struct seg {
  seg(int start, int end, double _slack) :
    viol(_slack){ends.first = start; ends.second = end;};
    
    int start() const {return ends.first;};
    int end() const {return ends.second;};

    bool operator <(const seg &val) const {
      return viol > val.viol;}
    
    std::pair<int, int> ends;
    double viol;
  };
  
  SupportGraph &G_s;
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;

  PSEPlp &m_lp;

  std::priority_queue<seg> pq; 
};

#endif
