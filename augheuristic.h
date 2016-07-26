#ifndef PSEP_AUG_HEURISTIC_H
#define PSEP_AUG_HEURISTIC_H

#include<iostream>
#include<vector>

#include "PSEP_util.h"
#include "lp.h"
#include "datagroups.h"

class PSEP_AugHeuristic {
 public:
 PSEP_AugHeuristic(PSEP_BestGroup &BestGroup, PSEP_LPGroup &LPGroup,
		   PSEP_SupportGroup &SupportGroup):
  best_tour_edges(BestGroup.best_tour_edges), m_lp(LPGroup.m_lp),
    support_indices(SupportGroup.support_indices),
    support_ecap(SupportGroup.support_ecap) {}
  
  int add_clamp();
  int clear_clamps();
  bool active(){ return !clamp_edges.empty();}

 private:
  std::vector<int> clamp_edges;

  std::vector<int> &best_tour_edges;
  PSEPlp &m_lp;
  std::vector<int> &support_indices;
  std::vector<double> &support_ecap;
};

#endif
