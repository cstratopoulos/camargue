#ifndef CMR_BBCONSTRAINTS_H
#define CMR_BBCONSTRAINTS_H

#include<array>
#include<memory>
#include<vector>

#include "lp.hpp"
#include "Graph.hpp"
#include "BButils.hpp"
#include "LPprune.hpp"
#include "LPcore.hpp"

namespace CMR {
  namespace BB {
    
    class Constraints {
    public:
    Constraints(CMR::BB::BranchPlan _Strategy,
		Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
		Data::LPGroup &LPGroup,
		Data::SupportGroup &SupportGroup,
		CMR::BB::RightBranch &_RB, CMR::BB::EdgeStatuses &_ES,
		CMR::LP::CutPrune &_LPPrune, CMR::LP::Core &_LPCore):
      ncount(GraphGroup.m_graph.node_count), Strategy(_Strategy),
	edges(GraphGroup.m_graph.edges),
	best_tour_edges(BestGroup.best_tour_edges),
	m_lp_edges(LPGroup.m_lp_edges), m_lp(LPGroup.m_lp),
	support_indices(SupportGroup.support_indices), RBranch(_RB),
	EdgeStats(_ES), LPPrune(_LPPrune), LPCore(_LPCore) {}

    public:
      int compute_branch_edge();
      int enforce(std::unique_ptr<CMR::BB::TreeNode> &v);
      int unenforce(std::unique_ptr<CMR::BB::TreeNode> &v);
      const int ncount;

    private:
      friend class CMR::BB::Visitor;

      CMR::BB::BranchPlan Strategy;

      bool naive_compatible(const int clamp, const int partner);
      
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
      CMRlp &m_lp;

      std::vector<int> &support_indices;

      std::vector<int> naive_branch_candidates;
      int naive_edge_partner;

      CMR::BB::RightBranch &RBranch;
      CMR::BB::EdgeStatuses &EdgeStats;

      CMR::LP::CutPrune &LPPrune;
      CMR::LP::Core &LPCore;
    };
  }
}

#endif
