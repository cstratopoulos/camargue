#ifndef PSEP_LP_PRUNE_H
#define PSEP_LP_PRUNE_H

#include<vector>

#include "lp.h"
#include "datagroups.h"
#include "PSEP_util.h"

namespace PSEP {
  namespace LP {
  class CutPrune {
  public:
  CutPrune(Data::GraphGroup &GraphGroup, Data::LPGroup &LPGroup) :
    m_lp(LPGroup.m_lp), node_count(GraphGroup.m_graph.node_count) {}

    int aggressive_prune();
    int prune_cuts(int &num_removed);
    int prune_with_skip(int &num_removed, IntPair skiprange,
			std::vector<int> &delset);
  private:    
    int dual_measure();
    
    PSEPlp &m_lp;
    int node_count;

    std::vector<double> dual_vars;
    std::vector<short int> ages;
    std::vector<int> delset;
  };
  }
}

#endif
