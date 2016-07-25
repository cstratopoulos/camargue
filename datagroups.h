#ifndef PSEP_DATAGROUP_H
#define PSEP_DATAGROUP_H

#include<vector>

#include "Graph.h"
#include "lp.h"

struct PSEP_GraphGroup {
  PSEP_GraphGroup(char *fname, CCdatagroup *dat);
  
  Graph m_graph;
  std::vector<int> island;
  std::vector<int> delta;
  std::vector<int> edge_marks;
};

struct PSEP_BestGroup {
  PSEP_BestGroup(const Graph &graph, CCdatagroup *dat);
  std::vector<int> best_tour_edges;
  std::vector<int> best_tour_nodes;
  std::vector<int> perm;

  double min_tour_value;
};

struct PSEP_SupportGroup {
  SupportGraph G_s;
  std::vector<int> support_indices;
  std::vector<int> support_elist;
  std::vector<double> support_ecap;
};

struct PSEP_LPGroup {
  PSEPlp m_lp;  
  std::vector<double> m_lp_edges;
};

#endif
