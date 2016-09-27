#ifndef PSEP_SEGCUTS_H
#define PSEP_SEGCUTS_H

#include<vector>

#include "lp.hpp"
#include "cuts.hpp"
#include "Graph.hpp"

namespace PSEP {
  template<> class Cut<seg> {
  public:
    Cut<seg>(std::vector<int> &_edge_marks, std::vector<int> &_best_tour_nodes,
	     SupportGraph &_G_s, CutQueue<HyperGraph> &segment_queue):
    edge_marks(_edge_marks), best_tour_nodes(_best_tour_nodes), G_s(_G_s),
      local_q(seg_q_max),
      external_q(segment_queue){}

    int cutcall();

  protected:
    int separate();
    int add_cuts();
    
  private:
    int build_hypergraph(const seg& seg_cut);
    
    std::vector<int> &edge_marks;
    const std::vector<int> &best_tour_nodes;
    
    SupportGraph &G_s;

    CutQueue<seg> local_q;
    CutQueue<HyperGraph> &external_q;
  };
}

#endif
