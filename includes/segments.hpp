#ifndef PSEP_SEGCUTS_H
#define PSEP_SEGCUTS_H

#include "cuts.hpp"
#include "Graph.hpp"

#include<vector>

namespace PSEP {

template<> class Cut<seg> {
public:
  Cut<seg>(std::vector<int> &_best_tour_nodes, std::vector<int> &_perm,
	   SupportGraph &_G_s, std::vector<int> &_support_elist,
	   std::vector<double> &_support_ecap,
	   CutQueue<HyperGraph> &segment_queue):
    best_tour_nodes(_best_tour_nodes), perm(_perm),
    G_s(_G_s), support_elist(_support_elist), support_ecap(_support_ecap),
    local_q(segment_queue.q_capacity),
    external_q(segment_queue){}

  int cutcall();

protected:
  int separate();
  int add_cuts();
    
private:
  int build_hypergraph(const seg& seg_cut);
  static int linsub_all_callback(double cut_val, int cut_start, int cut_end,
				 void *cb_data);
  static int linsub_callback(double vut_val, int cut_start, int cut_end,
			     void *cb_data);

  std::vector<int> &best_tour_nodes;
  std::vector<int> &perm;

  SupportGraph &G_s;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  CutQueue<seg> local_q;
  CutQueue<HyperGraph> &external_q;
};

}

#endif
