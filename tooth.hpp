#ifndef PSEP_TOOTH_HPP
#define PSEP_TOOTH_HPP

#include "Graph.hpp"

#include <vector>

namespace PSEP {

/*
 * Structure for storing the simple teeth which may comprise simple domino
 * parity inequalities. A simple tooth is a root node and a segment of a
 * best known tour, stored as best tour node indices. Hence the root 
 * node is best_tour_nodes[root], and the segment is
 * best_tour_nodes[body_start], best_tour_nodes[body_start + 1], ...
 *                                best_tour_nodes[body_end],
 * incrementing indices mod ncount. 
 * A SimpleTooth is meaningless without management by its CandidateTooth
 * class, which is responsible for all nontrivial operations on the tooth. 
 */

struct SimpleTooth {
  SimpleTooth(int _root, int _body_start, int _body_end) :
    root(_root), body_start(_body_start), body_end(_body_end) {}

  int root;
  int body_start;
  int body_end;

  //returns true iff root lies in the middle of the segment from body_start
  // to body_end
  bool sandwich() const;
};

class CandidateTeeth {
public:
  CandidateTeeth(std::vector<int> &_edge_marks,
		 std::vector<int> &_best_tour_nodes,
		 SupportGraph &_G_s) :
    edge_marks(_edge_marks), best_tour_nodes(_best_tour_nodes), G_s(_G_s) {}

  int body_size(const SimpleTooth &T);
  void complement(SimpleTooth &T);
  int body_subset(const SimpleTooth &T, const SimpleTooth &R, bool &result);

  void increment_slack(const int new_vx, double &lhs, int &rhs);
  
private:
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;
  SupportGraph &G_s;
};

}


#endif
