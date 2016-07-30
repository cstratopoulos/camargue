#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include "datagroups.h"
#include "segments.h"
#include "blossoms.h"
#include "cuts.h"

namespace PSEP{
  class CutControl {
  public:
  CutControl(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
		  PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup):
    //separate cut classes
    segments(GraphGroup.delta, GraphGroup.edge_marks, GraphGroup.m_graph.edges,
	     BestGroup.best_tour_nodes, LPGroup.m_lp, SupportGroup.G_s),
      blossoms(GraphGroup.delta, GraphGroup.edge_marks,
	       GraphGroup.m_graph.edges, BestGroup.best_tour_edges,
	       LPGroup.m_lp, LPGroup.m_lp_edges, SupportGroup.support_indices,
	       SupportGroup.support_elist, SupportGroup.support_ecap){}
	       

    PSEP::Cut<PSEP::seg> segments;
    PSEP::Cut<PSEP::blossom> blossoms;
  };
}


#endif
