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

    typedef TreeNode::NType NodeType;
    typedef TreeNode::NStat NodeStat;

    /***********************************************************************/

    class EdgeStatuses {
    public:
      EdgeStatuses(std::vector<double> &lower_bounds) {
	for(int i = 0; i < lower_bounds.size(); i++)
	  if(lower_bounds[i] == 1.0) FixedUp.insert(i);
      }

    private:  
      void add_branch_var(const std::unique_ptr<PSEP::BB::TreeNode> &v){
	if(v->type() == BB::NodeType::LEFT)
	  Left.insert(v->edge());
	else
	  Right.insert(v->edge());
      }

      void remove_branch_var(const std::unique_ptr<PSEP::BB::TreeNode> &v){
	if(v->type() == BB::NodeType::LEFT)
	  Left.erase(v->edge());
	else
	  Right.erase(v->edge());
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
  }
}

#endif
