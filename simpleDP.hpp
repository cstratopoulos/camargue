#ifndef PSEP_SIMPLEDP_H
#define PSEP_SIMPLEDP_H

#include "cuts.hpp"
#include "DPgraph.hpp"

namespace PSEP {

template<> class PSEP::Cut<PSEP::dominoparity> {
public:

  int cutcall();

protected:
  int separate();
  int add_cuts();

private:
  int build_hypergraph(const PSEP::dominoparity &dp_cut);

  std::vector<int> &delta;
  std::vector<int> &edge_marks;

  std::vector<int> &best_tour_nodes;
  std::vector<int> &perm;

  PSEP::SupportGraph &G_s;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  PSEP::CutQueue<PSEP::dominoparity> local_q;
}

}

#endif
