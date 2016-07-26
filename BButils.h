#ifndef PSEP_BBUTILS_H
#define PSEP_BBUTILS_H

#include<memory>
#include<vector>
#include<unordered_set>
#include<map>
#include<utility>

#include "lp.h"
#include "PSEP_util.h"

class PSEP_BBNode {
 public:
  enum class NType {
    ROOT, LEFT, RIGHT
  };

  enum class NStat {
    UNVISITED, VISITED, FATHOMED
  };

 PSEP_BBNode(NType _type, int _edge) :
  node_type(_type), node_stat(NStat::UNVISITED), branch_edge(_edge) {}

 PSEP_BBNode() : node_type(NType::ROOT), branch_edge(-1) {}

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

class PSEP_EdgeStatuses {
 public:
  PSEP_EdgeStatuses() {}
  
  void add_branch_var(const std::unique_ptr<PSEP_BBNode> &v){
    if(v->type() == PSEP_BBNode::NType::LEFT)
      Left.insert(v->edge());
    else
      Right.insert(v->edge());
  }
  void augmentation_merge(){
    Left.insert(Right.begin(), Right.end());
    Right.clear();
  }

 private:
  std::unordered_set<int> Left;
  std::unordered_set<int> Right;		 
};

class PSEP_RightBranch {
 public:
 PSEP_RightBranch() : constraint_range(IntPair(0,0)),
    first_right(-1) {}  
  bool active() const { return !edge_row_lookup.empty(); }
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
  int first_right_edge();

 private:
  friend class PSEP_BBConstraints;
  std::map<int, int> edge_row_lookup;
  IntPair constraint_range;
  int first_right;  
};

class PSEP_BBConstraints {
 public:
  PSEP_BBConstraints(PSEP_BestGroup &BestGroup, PSEP_LPGroup &LPGroup,
		     PSEP_RightBranch &_RB, PSEP_EdgeStatuses &_ES):
  best_tour_edges(BestGroup.best_tour_edges), m_lp_edges(LPGroup.m_lp_edges),
    m_lp(LPGroup.m_lp), RightBranch(_RB), EdgeStats(_ES) {}
  
  int compute_branch_edge();
  
  int add_first_right(const int edge);
  int update_right(const int edge);
  int remove_right(const int edge);
  
  int add_left(const int edge);
  int remove_left(const int edge);

 private:
  std::vector<int> &best_tour_edges;
  std::vector<double> &m_lp_edges;
  PSEPlp &m_lp;

  PSEP_RightBranch &RightBranch;
  PSEP_EdgeStatuses &EdgeStats;
};

#endif
