#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include<iostream>
#include<iomanip>

#include "datagroups.h"
#include "segments.h"
#include "blossoms.h"
#include "dominos.h"
#include "gencuts.h"
#include "safegmi.h"
#include "cuts.h"
#include "PSEP_util.h"

namespace PSEP{
  class CutControl {
  public:
  CutControl(Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
	     Data::LPGroup &LPGroup, Data::SupportGroup &SupportGroup):
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
      general_cuts(true, true, true, BestGroup.best_tour_edges, LPGroup.m_lp,
		   LPGroup.m_lp_edges, LPGroup.frac_colstat,
		   LPGroup.frac_rowstat),
      safe_gomory(BestGroup.best_tour_edges,
		  LPGroup.m_lp, LPGroup.m_lp_edges, LPGroup.frac_colstat,
		  LPGroup.frac_rowstat,
		  SupportGroup.support_indices),
      prefs(LPGroup.prefs),
      total_segtime(0), total_2mtime(0), total_dptime(0), total_gentime(0),
      total_segcalls(0), total_2mcalls(0), total_gencalls(0){}

    
  public:
    int primal_sep(const int augrounds, const LP::PivType stat);
    int general_sep(const double piv_val);
    int safe_gomory_sep();
    
    void profile(){
      std::cout << "   Total time during lightDP sep: " << std::setprecision(4)
		<< total_dptime << "s\n"
		<< "                     segment sep: "
		<< total_segtime << "s\n"
		<< "                     blossom sep: "
		<< total_2mtime << "s\n"
		<< "                         MIP sep: "
		<< total_gentime << "s\n" << std::setprecision(6);
    }

  private:
    PSEP::Cut<PSEP::seg> segments;
    PSEP::Cut<PSEP::blossom> blossoms;
    PSEP::Cut<PSEP::domino> dominos;

    PSEP::Cut<PSEP::general> general_cuts;

    PSEP::Cut<PSEP::safeGMI> safe_gomory;

    PSEP::LP::Prefs &prefs;

    double total_segtime, total_2mtime, total_dptime, total_gentime;
    int total_segcalls, total_2mcalls, total_gencalls;
  };
}


#endif
