#ifndef PSEP_CUTCALL_H
#define PSEP_CUTCALL_H

#include <vector>
#include <memory>

#include "lp.hpp"
#include "datagroups.hpp"
#include "segments.hpp"
#include "blossoms.hpp"
#include "fastblossoms.hpp"
#include "simpleDP.hpp"
#include "safegmi.hpp"
#include "tooth.hpp"
#include "cuts.hpp"
#include "PSEP_util.hpp"
#include "timer.hpp"


namespace PSEP{

class CutControl {
public:
  CutControl(Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
	     Data::LPGroup &LPGroup, Data::SupportGroup &SupportGroup,
	     const PSEP::Timer *purecut_timer_p):
    set_repo(BestGroup.best_tour_nodes, BestGroup.perm),
    translator(GraphGroup),
    segment_q(LPGroup.prefs.max_per_round),
    blossom_q(LPGroup.prefs.q_max_size),
    domino_q(10),
    segments(BestGroup.best_tour_nodes, BestGroup.perm,
	     SupportGroup.G_s, SupportGroup.support_elist,
	     SupportGroup.support_ecap,
	     segment_q),
    blossoms(GraphGroup.delta, GraphGroup.edge_marks,
	     GraphGroup.m_graph, BestGroup.best_tour_edges,
	     LPGroup.m_lp_edges, SupportGroup.support_indices,
	     SupportGroup.support_elist, SupportGroup.support_ecap,
	     blossom_q),
    fastblossoms(GraphGroup.delta, GraphGroup.edge_marks,
		 GraphGroup.m_graph, BestGroup.best_tour_edges,
		 LPGroup.m_lp_edges,
		 SupportGroup.support_indices,
		 SupportGroup.support_elist, SupportGroup.support_ecap,
		 blossom_q),
    safe_gomory(BestGroup.best_tour_edges,
		LPGroup.m_lp, LPGroup.m_lp_edges, LPGroup.frac_colstat,
		LPGroup.frac_rowstat,
		SupportGroup.support_indices,
		LPGroup.prefs.max_per_round),
    graph_data(GraphGroup),
    LP_data(LPGroup),
    supp_data(SupportGroup),
    best_data(BestGroup),
    segtime("Segment sep", purecut_timer_p),
    matchtime("2match sep", purecut_timer_p),
    dptime("Simple DP sep", purecut_timer_p),
    gmitime("Safe GMI sep", purecut_timer_p),
    total_segcalls(0), total_2mcalls(0), total_gencalls(0){}

    
public:
  int primal_sep(const int augrounds, const LP::PivType stat);
  int add_primal_cuts();
  
  int safe_gomory_sep();

  int in_subtour_poly(bool &result);
    
  void profile();

private:
  int q_has_viol(bool &result, CutQueue<HyperGraph> &pool_q);

  
  PSEP::SetBank set_repo;
  PSEP::CutTranslate translator;
  
  PSEP::CutQueue<PSEP::HyperGraph> segment_q;
  PSEP::CutQueue<PSEP::HyperGraph> blossom_q;

  PSEP::CutQueue<PSEP::dominoparity> domino_q;
  
  PSEP::Cut<PSEP::seg> segments;
  PSEP::Cut<PSEP::blossom> blossoms;
  PSEP::Cut<PSEP::fastblossom> fastblossoms;

  /** @todo: make all cut classes like this. */
  std::unique_ptr<PSEP::Cut<PSEP::dominoparity>> dominos;

  PSEP::Cut<PSEP::safeGMI> safe_gomory;

  PSEP::Data::GraphGroup &graph_data;
  PSEP::Data::LPGroup &LP_data;
  PSEP::Data::SupportGroup &supp_data;
  PSEP::Data::BestGroup &best_data;

  PSEP::Timer segtime, matchtime, dptime, gmitime;
  int total_segcalls, total_2mcalls, total_gencalls;
};

}


#endif
