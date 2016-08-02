#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include<iostream>
#include<iomanip>

#include "datagroups.h"
#include "segments.h"
#include "blossoms.h"
#include "dominos.h"
#include "cuts.h"
#include "PSEP_util.h"

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
	       SupportGroup.support_elist, SupportGroup.support_ecap),
      dominos(GraphGroup.edge_marks,
	      GraphGroup.m_graph.edge_lookup,
	      BestGroup.best_tour_nodes, BestGroup.perm, LPGroup.m_lp,
	      LPGroup.m_lp_edges, SupportGroup.G_s, SupportGroup.support_elist,
	      SupportGroup.support_ecap),
      prefs(LPGroup.prefs),
      total_segtime(0), total_2mtime(0), total_dptime(0),
      total_segcalls(0), total_2mcalls(0){}

    
  public:
    int primal_sep(const int augrounds, const PivType stat);
    
    void profile(){
      std::cout << "   Total time during lightDP sep: " << std::setprecision(4)
		<< total_dptime << "s\n";
      std::cout << "   Average time per segment call: "
		<< ((double) (total_segtime / total_segcalls)) << "\n";
      std::cout << "                     2match call: "
		<< ((double) (total_2mtime / total_2mcalls)) << "\n"
		<< std::setprecision(6);
    }

  private:
    PSEP::Cut<PSEP::seg> segments;
    PSEP::Cut<PSEP::blossom> blossoms;
    PSEP::Cut<PSEP::domino> dominos;

    PSEP_LP_Prefs &prefs;

    double total_segtime, total_2mtime, total_dptime;
    int total_segcalls, total_2mcalls;
  };
}


#endif
