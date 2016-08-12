#ifndef PSEP_SIMPLEDP_H
#define PSEP_SIMPLEDP_H

extern "C" {
#include <concorde/concorde.h>
}

#include "tooth.h"
#include "lp.h"
#include "cuts.h"

namespace PSEP {
  class SimpleDP {
  public:
  SimpleDP(std::vector<int> & _tour_nodes,
	   std::vector<int> & _perm, SupportGraph & _G,
	   std::vector<int> & _edge_marks, std::vector<int> & _elist,
	   std::vector<double>& _ecap, IntPairMap & _edge_lookup,
	   PSEPlp & _lp) :
    G_s(_G), best_tour_nodes(_tour_nodes), perm(_perm),
      support_elist(_elist), support_ecap(_ecap), edge_lookup(_edge_lookup),
      m_lp(_lp), candidates(_tour_nodes, _G, _edge_marks) {}
		
    int separate();
    void parse_domino(const int deltacount, const std::vector<int> &dom_delta,
		      std::vector<double> &agg_coeffs, double *rhs_p);
    void print_cutgraph(const int ncount, const int webcount);

  private:
    void build_light_cuttree();
    void add_web_edges();
    int call_CC_gomoryhu();
    int in_subtour_poly(bool *result_p);
  
    friend class PSEP::Cut<PSEP::domino>;
    SupportGraph &G_s;

    std::vector<int> &best_tour_nodes;
    std::vector<int> &perm;

    std::vector<int> &support_elist;
    std::vector<double> &support_ecap;

    IntPairMap &edge_lookup;

    PSEPlp &m_lp;

    std::vector<std::shared_ptr<CandTooth::SimpleTooth>> light_nodes;
    std::vector<int> cut_elist;
    std::vector<double> cut_ecap;
    std::vector<int> cut_marks;
    std::vector<bool> node_marks;

    std::vector<int> cut_nodes;
  
    CandTooth candidates;
  };
}

#endif
