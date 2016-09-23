#ifndef PSEP_LP_FIXING_H
#define PSEP_LP_FIXING_H

#include<vector>

#include "lp.h"
#include "datagroups.h"
#include "Graph.h"
#include "PSEP_util.h"

namespace PSEP {
  namespace LP {
class EdgeFix {
 public:
 EdgeFix(PSEP::Data::BestGroup &BestGroup,
	    PSEP::Data::GraphGroup &GraphGroup,
	    PSEP::Data::LPGroup &LPGroup) :
  m_min_tour_value(BestGroup.min_tour_value),
    best_tour_edges(BestGroup.best_tour_edges),
    m_graph(GraphGroup.m_graph), edges(GraphGroup.m_graph.edges),
    edge_lookup(GraphGroup.m_graph.edge_lookup), delta(GraphGroup.delta),
    m_lp(LPGroup.m_lp), m_lp_edges(LPGroup.m_lp_edges)
    {}

  enum FixStats {
    LEAVE = 0, DELETE = 1, FIXED = 2
  };

  int price(int *clamptotal, int *deltotal);
  int fixup();
  int delete_cols();
  void delete_edges();

  int redcost_fixing();
  

 private:
  double &m_min_tour_value;
  std::vector<int> &best_tour_edges;
  
  Graph &m_graph;
  std::vector<Edge> &edges;
  IntPairMap &edge_lookup;
  std::vector<int> &delta;

  PSEPlp &m_lp;
  std::vector<double> &m_lp_edges;

  std::vector<int> edge_delset;

};
  }
}

#endif
