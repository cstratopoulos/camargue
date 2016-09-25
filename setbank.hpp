#ifndef PSEP_SETBANK_HPP
#define PSEP_SETBANK_HPP

#include<unordered_set>
#include<vector>
#include<iostream>

#include<boost/functional/hash.hpp>

#include "PSEP_util.hpp"

namespace PSEP{

void hashing_test();
void interactive_test();

struct IntervalSet {
  IntervalSet(std::vector<int> &nodelist,
	      const std::vector<int> &best_tour_nodes,
	      const std::vector<int> &perm);
  IntervalSet(std::vector<IntPair> _interval_list) :
    interval_list(_interval_list), use_count(0) {}
  
  std::vector<IntPair> interval_list;
  int use_count;

  void add_use() const { const_cast<IntervalSet *>(this)->use_count += 1; }
  void del_use() const { const_cast<IntervalSet *>(this)->use_count -= 1; }
  
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
typedef SetHash::iterator SetRef;

struct HyperGraph {
  std::vector<SetRef> set_refs;
  int rhs;
  int row_number;
};

class SetBank {
  SetBank(std::vector<int> &best_tour_nodes, std::vector<int> &_perm) :
    tour_nodes(best_tour_nodes), perm(_perm) {}

  void get_hypergraph(HyperGraph &newcut,
		      std::vector<std::vector<int>> &hypergraph_edges);
  void delete_hypergraph(HyperGraph &oldcut);  
  
private:
  SetRef add_or_increment(std::vector<int> nodeset);
  void remove_or_decrement(std::vector<int> nodeset);
  
  SetHash set_bank;
  
  std::vector<int> tour_nodes;
  std::vector<int> perm;
};


}

#endif
