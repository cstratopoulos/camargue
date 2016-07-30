#ifndef PSEP_BLOSSOMS_H
#define PSEP_BLOSSOMS_H

#include<vector>

#include "cuts.h"
#include "Graph.h"
#include "lp.h"

namespace PSEP {
  struct blossom {
  blossom(std::vector<int> _handle, int _cut_edge, double _val) :
    handle(_handle), cut_edge(_cut_edge), cut_val(_val){}

    bool operator< (const blossom &val) const {
      return cut_val < val.cut_val;
    }

    std::vector<int> handle;
    int cut_edge;
    double cut_val;
  };
}

namespace PSEP{
  template<> class Cut<blossom> {
  public:
    Cut<blossom>(std::vector<int> &_delta, std::vector<int> &_edge_marks,
		 std::vector<Edge> &_edges, std::vector<int> &_best_tour_edges,
		 PSEPlp &_m_lp, std::vector<double> &_m_lp_edges,
		 std::vector<int> &_support_indices,
		 std::vector<int> &_support_elist,
		 std::vector<double> &_support_ecap):
    deltacount(0), delta(_delta), edge_marks(_edge_marks), edges(_edges),
      best_tour_edges(_best_tour_edges), m_lp(_m_lp),
      m_lp_edges(_m_lp_edges),
      support_indices(_support_indices),
      support_elist(_support_elist), support_ecap(_support_ecap) {}

    int cutcall();

  private:
    int separate();
    int parse_coeffs();
    int add_cut();

    int deltacount;
    std::vector<int> &delta;
    std::vector<int> &edge_marks;
    std::vector<Edge> &edges;
    std::vector<int> &best_tour_edges;
    PSEPlp &m_lp;
    std::vector<double> &m_lp_edges;
    std::vector<int> &support_indices;
    std::vector<int> &support_elist;
    std::vector<double> &support_ecap;

    std::vector<double> cut_ecap;

    //    G_Utils gu;

    std::unique_ptr<blossom> best;
  };
}

#endif
