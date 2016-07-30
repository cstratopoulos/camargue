#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include "datagroups.h"
#include "cuts.h"
#include "segments2.h"

namespace PSEP{
  class CutControl {
  public:
  CutControl(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
		  PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup):
    //separate cut classes
    segments(GraphGroup.delta, GraphGroup.edge_marks,
	     GraphGroup.m_graph.edges,
	     BestGroup.best_tour_nodes,
	     LPGroup.m_lp, SupportGroup.G_s){}
	       

    PSEP::Cut<PSEP::seg> segments;
  };
}


#endif
