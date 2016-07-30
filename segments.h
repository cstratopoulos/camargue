#ifndef PSEP_SEGCUTS_H
#define PSEP_SEGCUTS_H

#include<vector>

#include "cuts.h"
#include "Graph.h"
#include "lp.h"

namespace PSEP {
  struct seg {
  seg(int _start, int _end, double _slack) :
    start(_start), end(_end), viol(_slack) {}
    int start;
    int end;
    double viol;

    bool operator <(const seg &val) const {
      return viol > val.viol;
    }
  };
}

namespace PSEP {
  template<> class Cut<seg> {
  public:
    Cut<seg>(std::vector<int> &_delta, std::vector<int> &_edge_marks,
	     std::vector<Edge> &_edges, std::vector<int> &_best_tour_nodes,
	     PSEPlp &_m_lp, SupportGraph &_G_s):
    deltacount(0), delta(_delta), edge_marks(_edge_marks), edges(_edges),
      best_tour_nodes(_best_tour_nodes), m_lp(_m_lp), G_s(_G_s) {}

    int cutcall();

  private:
    int separate();
    int parse_coeffs();
    int add_cut();

    int deltacount;
    std::vector<int> &delta;
    std::vector<int> &edge_marks;
    /*const */std::vector<Edge> &edges;
    const std::vector<int> &best_tour_nodes;
    PSEPlp &m_lp;

    SupportGraph &G_s;

    G_Utils gu;

    std::unique_ptr<seg> best;
  };
}

#endif
