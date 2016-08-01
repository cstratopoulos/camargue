#ifndef PSEP_BBUTILS_H
#define PSEP_BBUTILS_H

#include<array>
#include<memory>
#include<vector>
#include<unordered_set>
#include<map>
#include<utility>

#include "datagroups.h"
#include "lp.h"
#include "PSEP_util.h"

namespace PSEP {
  namespace BB {
    class TreeNode;
    class EdgeStatuses;
    class RightBranch;
    class Constraints;
    
    class TreeNode {
    public:
      enum class NType {
	ROOT, LEFT, RIGHT
	  };

      enum class NStat {
	UNVISITED, VISITED, FATHOMED
	  };

    TreeNode(NType _type, int _edge) :
      node_type(_type), node_stat(NStat::UNVISITED), branch_edge(_edge) {}

    TreeNode() : node_type(NType::ROOT), branch_edge(-1) {}

      NType type() const { return node_type; }
      NStat status() const { return node_stat; }

      int edge() const {
	if(node_type != NType::ROOT)
	  return branch_edge;
	else
	  return -1;
      }

    private:
      NType node_type;
      NStat node_stat;
      int branch_edge;
    };

    /***********************************************************************/

    class EdgeStatuses {
    public:
      EdgeStatuses(std::vector<double> &lower_bounds) {
	for(int i = 0; i < lower_bounds.size(); i++)
	  if(lower_bounds[i] == 1.0) FixedUp.insert(i);
      }

    private:  
      int add_branch_var(const std::unique_ptr<PSEP::BB::TreeNode> &v){
	if(v->type() == PSEP::BB::TreeNode::NType::ROOT){
	  std::cerr << "EdgeStatuses::add_branch_var tried to add root\n";
	  return 1;
	}

	if(v->status() != PSEP::BB::TreeNode::NStat::UNVISITED){
	  std::cerr << "EdgeStatuses::add_branch_var tried to add visited "
		    << "node\n"; return 1;
	}
	
	if(v->type() == PSEP::BB::TreeNode::NType::LEFT)
	  Left.insert(v->edge());
	else
	  Right.insert(v->edge());
	
	return 0;
      }

      int remove_branch_var(const std::unique_ptr<PSEP::BB::TreeNode> &v){
	if(v->type() == PSEP::BB::TreeNode::NType::ROOT){
	  std::cerr << "EdgeStatuses::remove_branch_var tried to add root\n";
	  return 1;
	}

	if(v->status() == PSEP::BB::TreeNode::NStat::UNVISITED){
	  std::cerr << "EdgeStatuses::remove_branch_var tried to remove "
	       << "unvisited node\n"; return 1;
	}
	
	if(v->type() == PSEP::BB::TreeNode::NType::LEFT)
	  Left.erase(v->edge());
	else
	  Right.erase(v->edge());

	return 0;
      }
      
      void augmentation_merge(){
	Left.insert(Right.begin(), Right.end());
	Right.clear();
      }
      
      friend class PSEP::BB::RightBranch;
      friend class PSEP::BB::Constraints;
      std::unordered_set<int> Left;
      std::unordered_set<int> Right;
      std::unordered_set<int> FixedUp;
    };
    
    /***********************************************************************/

    class RightBranch {
    public:
    RightBranch() : constraint_range(IntPair(0,0)),
	first_right(-1) {}  
      bool active() const { return first_right != -1; }
      int row_lookup(const int edge) const {
	std::map<int, int>::const_iterator
	  it = edge_row_lookup.find(edge);
	if(it != edge_row_lookup.end())
	  return it->second;
	else
	  return -1;
      }
  
      void get_range(int *first, int *last) const {
	*first = constraint_range.first; *last = constraint_range.second;
      }
      
      void update_range(const std::vector<int> &delset){
	for(std::map<int, int>::iterator it = edge_row_lookup.begin();
	    it != edge_row_lookup.end(); it++)
	  it->second = delset[it->second];
      }
      
      void add_edge_row(const int edge, const int rownum){
	edge_row_lookup[edge] = rownum;
      }
      
      int first_right_edge() const {return first_right;} 

    private:
      friend class PSEP::BB::Constraints;
      std::map<int, int> edge_row_lookup;
      IntPair constraint_range;
      int first_right;  
    };

    /***********************************************************************/

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

      int add_right_branch(const int edge);
      int add_first_right_rows(const int edge);
      int explore_right(const int edge);
      int update_right(const int edge);
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
