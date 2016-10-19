#ifndef PSEP_TOOTH_HPP
#define PSEP_TOOTH_HPP

#include "Graph.hpp"

#include <memory>
#include <vector>

#define PSEP_TOOTH_UNIQ


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
    root(_root), body_start(_body_start), body_end(_body_end),
    slack(INFINITY) {}

  #ifdef PSEP_TOOTH_UNIQ
  typedef std::unique_ptr<SimpleTooth> Ptr;
  #else
  typedef SimpleTooth* Ptr;
  #endif

  int root;
  int body_start;
  int body_end;

  double slack;

  //returns true iff root lies in the middle of the segment from body_start
  // to body_end
  bool sandwich() const;
};

/*
 * This class generates and manages a collection of candidate teeth to be 
 * included in a simple domino parity inequality. It is also responsible for
 * nontrivial operations on teeth which require knowledge of the best tour
 * and perm and underlying graph used to define them. 
 */

class CandidateTeeth {
public:
  CandidateTeeth(std::vector<int> &_edge_marks,
		 std::vector<int> &_best_tour_nodes,
		 SupportGraph &_G_s) :
    light_teeth(std::vector<std::vector<SimpleTooth::Ptr>>(_G_s.node_count)),
    //TODO: this may throw
    edge_marks(_edge_marks), best_tour_nodes(_best_tour_nodes), G_s(_G_s) {}

  int body_size(const SimpleTooth &T);

  //if T has root i, body S, T complement has root i, body V\(Su{i})
  void complement(SimpleTooth &T);

  /*
   * For teeth T, R with the same root, returns whether the body of T is a
   * subset of the body of R, storing the result in result.
   * TODO: generalize to arbitrary teeth
   */
  int body_subset(const SimpleTooth &T, const SimpleTooth &R, bool &result);


  /*
   * This function is used to incrementally compute the slack on a collection
   * of SimpleTooth structs with the same root. 
   * Given a root i, we first construct the tooth T with trivial body equal
   * to a single edge from the current best tour.
   * lhs is initialized to zero, rhs is initialized to -1
   * The assumption is that T will be `incremented' to a tooth with the same
   * root, and body equal to the next node, new_vx, in the tour. Note that
   * new_vx is some actual node index, NOT best_tour_nodes[new_vx].
   * This routine will mark a single node in edge_marks, and it is the 
   * responsibility of the calling routine to zero it out appropriately
   * before and after.
   */
  void increment_slack(SimpleTooth &T, const int new_vx,
		       double &lhs, int &rhs);

  /*
   * Build a collection of light teeth: teeth for which the slack of the 
   * associated simple tooth inequality is less than 1/2
   */
  int get_light_teeth();

  void print_tooth(const SimpleTooth &T);
  void print_collection();
  
  std::vector<std::vector<SimpleTooth::Ptr>> light_teeth;

  
private:
  void clear_collection();
  int get_adjacent_teeth(const int root);
  int get_distant_teeth(const int root);
  
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;
  SupportGraph &G_s;
};

}


#endif
