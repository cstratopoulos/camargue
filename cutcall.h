#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include "datagroups.h"
#include "cuts.h"
#include "blossom.h"
#include "segments2.h"
#include "simpleDP.h"

class PSEP_Cutcall {
 public:
  PSEP_Cutcall(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
	       PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup):
  segments(GraphGroup, BestGroup, SupportGroup, LPGroup),
  //GraphGroup
  edges(GraphGroup.m_graph.edges),
    delta(GraphGroup.delta),
    edge_marks(GraphGroup.edge_marks),
    //BestGroup
    best_tour_nodes(BestGroup.best_tour_nodes),
    //SupportGroup
    support_elist(SupportGroup.support_elist),
    support_ecap(SupportGroup.support_ecap),
    //LPGroup
    m_lp_edges(LPGroup.m_lp_edges),
    //separate cut classes
    blossoms(BestGroup.best_tour_edges, LPGroup.m_lp_edges,
	     SupportGroup.support_indices, SupportGroup.support_elist,
	     SupportGroup.support_ecap, LPGroup.m_lp),
    dominos(BestGroup.best_tour_nodes, BestGroup.perm, SupportGroup.G_s,
	    GraphGroup.edge_marks, SupportGroup.support_elist,
	    SupportGroup.support_ecap, GraphGroup.m_graph.edge_lookup,
	    LPGroup.m_lp){}
	       

  PSEP::Cuts<PSEP::seg> segments;
  int blossom(const int max_cutcount, int *num_added_p);
  int simpleDP(const int max_cutcount, int *num_added_p, int *num_bad_p);

  int in_subtour_poly(bool *result_p);

 private:
  std::vector<Edge> &edges;
  std::vector<int> &delta;
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;

  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  std::vector<double> &m_lp_edges;
   
  PSEP_2match blossoms;
  PSEP_SimpleDP dominos;
  
};


#endif
