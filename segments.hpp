#ifndef PSEP_SEGCUTS_H
#define PSEP_SEGCUTS_H

#include<vector>

#include "lp.hpp"
#include "cuts.hpp"
#include "Graph.hpp"

namespace PSEP {
  template<> class Cut<seg> {
  public:
    Cut<seg>(std::vector<int> &_delta, std::vector<int> &_edge_marks,
	     std::vector<Edge> &_edges, std::vector<int> &_best_tour_nodes,
	     PSEPlp &_m_lp, SupportGraph &_G_s,
	     CutQueue<HyperGraph> &segment_queue):
    deltacount(0), delta(_delta), edge_marks(_edge_marks), edges(_edges),
      best_tour_nodes(_best_tour_nodes), m_lp(_m_lp), G_s(_G_s),
      local_q(seg_q_max),
      external_q(segment_queue){}

    int cutcall();

  protected:
    int separate();
    int add_cut();
    
  private:
    int parse_coeffs();
    
    int deltacount;
    std::vector<int> &delta;
    std::vector<int> &edge_marks;
    /*const */std::vector<Edge> &edges;
    const std::vector<int> &best_tour_nodes;
    PSEPlp &m_lp;

    SupportGraph &G_s;

    std::unique_ptr<seg> best;

    CutQueue<seg> local_q;
    CutQueue<HyperGraph> &external_q;
  };
}

#endif
