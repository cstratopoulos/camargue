#ifndef PSEP_PURE_CUT_H
#define PSEP_PURE_CUT_H

#include <iostream>
#include <vector>

#include "datagroups.h"
#include "Graph.h"
#include "lp.h"
#include "PSEP_util.h"
#include "cutcall.h"
#include "LPcore.h"
#include "printer.h"

class PSEP_PureCut {
 public:
  PSEP_PureCut(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
	       PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup):
  /* //GraphGroup */
  /* m_graph(GraphGroup.m_graph), island(GraphGroup.island), */
  /*   delta(GraphGroup.delta), edge_marks(GraphGroup.edge_marks), */
  /*   //BestGroup */
  /*   best_tour_edges(BestGroup.best_tour_edges), */
  /*   best_tour_nodes(BestGroup.best_tour_nodes), perm(BestGroup.perm), */
  /*   m_min_tour_value(BestGroup.min_tour_value), */
  /*   //LPgroup */
  /*   m_lp(LPGroup.m_lp), */
  /*   m_lp_edges(LPGroup.m_lp_edges), */
  /*   //SupportGroup */
  /*   G_s(SupportGroup.G_s), support_indices(SupportGroup.support_indices), */
  /*   support_elist(SupportGroup.support_elist), */
  /*   support_ecap(SupportGroup.support_ecap), */
    print(BestGroup.best_tour_nodes, BestGroup.best_tour_edges,
	  LPGroup.m_lp_edges, GraphGroup.m_graph.edges),
    cutcall(GraphGroup, BestGroup, LPGroup, SupportGroup),
    LPcore(LPGroup, GraphGroup, SupportGroup, BestGroup){}

  int solve();

 private:
  /* Graph &m_graph; */
  /* std::vector<int> &island; */
  /* std::vector<int> &delta; */
  /* std::vector<int> &edge_marks; */

  /* std::vector<int> &best_tour_edges; */
  /* std::vector<int> &best_tour_nodes; */
  /* //perm[best_tour_nodes[i]] == i, perm[j] returns index of j in best tour */
  /* std::vector<int> &perm; */
  /* double m_min_tour_value; */
  
  /* PSEPlp &m_lp;   */
  /* std::vector<double> &m_lp_edges; */
  
  /* SupportGraph &G_s; */
  /* std::vector<int> &support_indices; */
  /* std::vector<int> &support_elist; */
  /* std::vector<double> &support_ecap; */

  
  /* G_Utils gu; */
  PSEP_Printer print;

  PSEP_Cutcall cutcall;
  PSEP_LP_Core LPcore;
};

#endif
