#ifndef CMR_SETBANK_HPP
#define CMR_SETBANK_HPP

#include<unordered_set>
#include<vector>
#include<iostream>

#include<boost/functional/hash.hpp>

#include "util.hpp"

namespace CMR{

void interactive_test();

struct IntervalSet {
  IntervalSet(std::vector<int> &nodelist,
	      const std::vector<int> &best_tour_nodes,
	      const std::vector<int> &perm);
  IntervalSet(std::vector<IntPair> _interval_list) :
    interval_list(_interval_list), use_count(0) {}
  
  std::vector<IntPair> interval_list;
  int use_count;

  void add_use() const { const_cast<IntervalSet*>(this)->use_count += 1; }
  void del_use() const { const_cast<IntervalSet*>(this)->use_count -= 1; }
  
  explicit operator bool() const { return !(interval_list.empty()); }
  bool operator==(const IntervalSet &rhs) const {
    return interval_list == rhs.interval_list;
  }
  
  void print(){
    for(const IntPair &p : interval_list)
      std::cout << p.first << ", " << p.second << "\n";
  }
};

//this is needed so boost::hash can hash the struct
std::size_t hash_value(const IntervalSet &S);

typedef std::unordered_set<IntervalSet, boost::hash<IntervalSet>> SetHash;

//forward declaration for static unique pointer
class SetBank;
  
struct HyperGraph {
  enum class CutType {
    Segment, Blossom
  };
  
  HyperGraph(std::vector<std::vector<int>> &node_sets,
	     const CutType _cut_type);
  HyperGraph(const IntPair &segment_ends);
  static bool same_tour(const std::vector<int> &best_tour_nodes);

  void delete_refs();
  CutType cut_type;
  std::vector<IntervalSet*> set_refs;

  static SetBank *source_setbank;
};

class SetBank {
public:
  SetBank(std::vector<int> &best_tour_nodes, std::vector<int> &_perm);
  int extract_nodelist(const IntervalSet &interval_set,
		       std::vector<int> &result_nodelist);
  
private:
  friend HyperGraph::HyperGraph(std::vector<std::vector<int>> &node_sets,
				const CutType _cut_type);
  friend HyperGraph::HyperGraph(const IntPair &segment_ends);
  friend void HyperGraph::delete_refs();
  friend bool HyperGraph::same_tour(const std::vector<int> &best_tour_nodes);
  IntervalSet *add_or_increment(IntervalSet &newset);
  void del_or_decrement(IntervalSet &oldset);
  
  SetHash set_bank;
  
  std::vector<int> tour_nodes;
  std::vector<int> perm;
};


}

#endif
