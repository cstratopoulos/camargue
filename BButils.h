#ifndef PSEP_BBUTILS_H
#define PSEP_BBUTILS_H

#include<array>
#include<iostream>
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
    class Visitor;
    
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
      friend class PSEP::BB::Visitor;
      
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
    RightBranch() : constraint_range(IntPair(-1,-1)),
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

      const IntPair &skiprange = constraint_range;
  
      
      int update_range(const std::vector<int> &delset){
	for(std::map<int, int>::iterator it = edge_row_lookup.begin();
	    it != edge_row_lookup.end(); it++){
	  if(delset[it->second] == -1){
	    std::cerr << "DELETED A RIGHT BRANCH ROW!\n"; return 1;
	  }
	  it->second = delset[it->second];

	}

	constraint_range.first = delset[constraint_range.first];
	constraint_range.second = delset[constraint_range.second];
	if(constraint_range.first == -1 || constraint_range.second == -1){
	  std::cerr << "BEGIN OR END OF RIGHT CONSTRAINT RANGE WAS DELETED!!\n";
	  return 1;
	}

	return 0;
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
