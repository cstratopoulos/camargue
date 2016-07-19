#ifndef PSEP_SIMPLEDP_H
#define PSEP_SIMPLEDP_H

extern "C" {
#include "../programs/concorde/concorde.h"
}

#include "tooth.h"
#include "lp.h"

class PSEP_SimpleDP {
 public:
  PSEP_SimpleDP(std::vector<int> & _tour_nodes,
		std::vector<int> & _perm, SupportGraph & _G,
		std::vector<int> & _edge_marks, std::vector<int> & _elist,
		std::vector<double>& _ecap, IntPairMap & _edge_lookup,
		PSEPlp & _lp) :
  G_s(_G), best_tour_nodes(_tour_nodes), perm(_perm),
    support_elist(_elist), support_ecap(_ecap), edge_lookup(_edge_lookup),
    m_lp(_lp), candidates(_tour_nodes, _G, _edge_marks) {}
		
  int separate(const int max_cutcount);
  void build_light_cuttree();
  void add_web_edges();
  int call_CC_gomoryhu(const int max_cutcount);
  void parse_domino(const int deltacount, const std::vector<int> &dom_delta,
		    std::vector<double> &agg_coeffs, double *rhs_p);
  int add_cut(const std::vector<double> &agg_coeffs, const double RHS);

 private:
  friend class PSEP_Cutcall;
  SupportGraph &G_s;

  std::vector<int> &best_tour_nodes;
  std::vector<int> &perm;

  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  IntPairMap &edge_lookup;

  PSEPlp &m_lp;

  std::vector<std::shared_ptr<PSEP_CandTooth::SimpleTooth> > light_nodes;
  std::vector<int> cut_elist;
  std::vector<double> cut_ecap;
  std::vector<int> cut_marks;
  std::vector<bool> node_marks;

  std::vector<std::vector<int> > toothlists;
  
  PSEP_CandTooth candidates;
};

#endif
