#ifndef PSEP_BBCONSTRAINTS_H
#define PSEP_BBCONSTRAINTS_H

#include<array>
#include<memory>
#include<vector>

#include "lp.h"
#include "Graph.h"
#include "BButils.h"
#include "LPprune.h"
#include "LPcore.h"

namespace PSEP {
  namespace BB {
    
    class Constraints {
    public:
    Constraints(PSEP::BB::BranchPlan _Strategy,
		Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
		Data::LPGroup &LPGroup,
		Data::SupportGroup &SupportGroup,
		PSEP::BB::RightBranch &_RB, PSEP::BB::EdgeStatuses &_ES,
		PSEP::LP::CutPrune &_LPPrune, PSEP::LP::Core &_LPCore):
      ncount(GraphGroup.m_graph.node_count), Strategy(_Strategy),
	edges(GraphGroup.m_graph.edges),
	best_tour_edges(BestGroup.best_tour_edges),
	m_lp_edges(LPGroup.m_lp_edges), m_lp(LPGroup.m_lp),
	support_indices(SupportGroup.support_indices), RBranch(_RB),
	EdgeStats(_ES), LPPrune(_LPPrune), LPCore(_LPCore) {}

    public:
      int compute_branch_edge();
      int enforce(std::unique_ptr<PSEP::BB::TreeNode> &v);
      int unenforce(std::unique_ptr<PSEP::BB::TreeNode> &v);
      const int ncount;

    private:
      friend class PSEP::BB::Visitor;

      PSEP::BB::BranchPlan Strategy;
      
      int add_left_clamp(const int edge);
      int remove_left_clamp(const int edge);

      int compute_partner_edge(const int clamp);
      void compute_right_row(const int clamp, const int partner,
			     std::array<double, 2> &rmatval, double &RHS);

      int compute_right_update(const int clamp, const int partner,
			     std::array<double, 2> &rmatval, double &RHS,
			     const std::vector<double> &new_tour);

      int add_right_branch(const int edge);
      
      int add_main_right_rows(const int edge);
      int add_naive_right_row(const int edge);
      
      int explore_right(const int edge);
      
      int update_right_rows();
      int remove_right(const int edge);

      int prune();
  
      std::vector<Edge> &edges;
  
      std::vector<int> &best_tour_edges;
  
      std::vector<double> &m_lp_edges;
      PSEPlp &m_lp;

      std::vector<int> &support_indices;

      PSEP::BB::RightBranch &RBranch;
      PSEP::BB::EdgeStatuses &EdgeStats;

      PSEP::LP::CutPrune &LPPrune;
      PSEP::LP::Core &LPCore;
    };
  }
}

#endif
