#ifndef PSEP_BBCONSTRAINTS_H
#define PSEP_BBCONSTRAINTS_H

#include<array>
#include<memory>
#include<vector>

#include "lp.h"
#include "Graph.h"
#include "BButils.h"

namespace PSEP {
  namespace BB {
    class Constraints {
    public:
    Constraints(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
		PSEP_LPGroup &LPGroup,
		PSEP_SupportGroup &SupportGroup,
		PSEP::BB::RightBranch &_RB, PSEP::BB::EdgeStatuses &_ES):
      edges(GraphGroup.m_graph.edges),
	best_tour_edges(BestGroup.best_tour_edges),
	m_lp_edges(LPGroup.m_lp_edges), m_lp(LPGroup.m_lp),
	support_indices(SupportGroup.support_indices), RBranch(_RB),
	EdgeStats(_ES) {}

    public:
      int compute_branch_edge();
      int enforce(std::unique_ptr<PSEP::BB::TreeNode> &v);
      int unenforce(std::unique_ptr<PSEP::BB::TreeNode> &v);

    private:  
      int add_left_clamp(const int edge);
      int remove_left_clamp(const int edge);

      void compute_right_row(const int clamp, const int partner,
			     std::array<double, 2> &rmatval, double &RHS){
	double clamp_best = best_tour_edges[clamp],
	       partner_best = best_tour_edges[partner];
	RHS = clamp_best - partner_best;
	rmatval = {2 * clamp_best - 1, 1 - 2 * partner_best};
      }

      int compute_right_update(const int clamp, const int partner,
			     std::array<double, 2> &rmatval, double &RHS,
			     const std::vector<double> &new_tour);

      int add_right_branch(const int edge);
      int add_first_right_rows(const int edge);
      int explore_right(const int edge);
      int update_right_rows();
      int remove_right(const int edge);
  
      std::vector<Edge> &edges;
  
      std::vector<int> &best_tour_edges;
  
      std::vector<double> &m_lp_edges;
      PSEPlp &m_lp;

      std::vector<int> &support_indices;

      PSEP::BB::RightBranch &RBranch;
      PSEP::BB::EdgeStatuses &EdgeStats;
    };
  }
}

#endif
