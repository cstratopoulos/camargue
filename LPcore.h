#ifndef PSEP_LPWORK_H
#define PSEP_LPWORK_H

#include<iostream>

#include "lp.h"
#include "Graph.h"
#include "datagroups.h"
#include "LPPrune.h"

namespace PSEP {
  class PureCut;
  class CutControl;
}

class PSEP_LP_Core {
 public:
  PSEP_LP_Core(PSEP_LPGroup &LPGroup, PSEP_GraphGroup &GraphGroup,
	       PSEP_SupportGroup &SupportGroup,
	       PSEP_BestGroup &BestGroup,
	       PSEP::LPPrune &_LPPrune) :
  LPPrune(_LPPrune),
  m_lp(LPGroup.m_lp), m_graph(GraphGroup.m_graph), prefs(LPGroup.prefs),
    old_colstat(LPGroup.old_colstat), old_rowstat(LPGroup.old_rowstat),
    m_lp_edges(LPGroup.m_lp_edges), G_s(SupportGroup.G_s),
    support_indices(SupportGroup.support_indices),
    support_elist(SupportGroup.support_elist),
    support_ecap(SupportGroup.support_ecap),
    best_tour_edges(BestGroup.best_tour_edges),
    best_tour_nodes(BestGroup.best_tour_nodes), perm(BestGroup.perm),
    m_min_tour_value(BestGroup.min_tour_value),
    island(GraphGroup.island), delta(GraphGroup.delta),
    edge_marks(GraphGroup.edge_marks){
    basis_init();
    if(prefs.switching_choice == LP::PRICING::SWITCHING::START){
      std::cout << "Immediate: ";
      change_pricing();
    }
  }
    
  bool is_dual_feas();
  bool is_integral();

  int single_pivot();
  int nondegenerate_pivot();
  int pivot_back();
  int primal_opt();

  int rebuild_basis(bool prune);
  int rebuild_basis(int &numremoved, IntPair skiprange,
		    std::vector<int> &delset);
  int basis_init();
  
  int pivot_until_change(PSEP::PivType &pivot_status);


  double get_obj_val();
  double set_edges();
  double set_support_graph();

  int update_best_tour();
  bool test_new_tour();

  int numrows(){ return PSEPlp_numrows(&m_lp); }
  int numcols(){ return PSEPlp_numcols(&m_lp); }

  void change_pricing();

 private:
  PSEP::LPPrune &LPPrune;
  PSEPlp &m_lp;
  Graph &m_graph;

  PSEP_LP_Prefs prefs;

  friend class TSP_Solver;
  friend class PSEP::PureCut;
  friend class PSEP::CutControl;
  std::vector<int> &old_colstat;
  std::vector<int> &old_rowstat;

  std::vector<double> &m_lp_edges;

  SupportGraph &G_s;
  
  std::vector<int> &support_indices;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;


  std::vector<int> &best_tour_edges;
  std::vector<int> &best_tour_nodes;
  std::vector<int> &perm;

  double &m_min_tour_value;

  std::vector<int> &island;
  std::vector<int> &delta;
  std::vector<int> &edge_marks;
};


#endif
