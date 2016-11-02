/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                BUILDING COLLECTIONS OF CANDIDATE TEETH
 *
 * This file contains structures and classes used to build a collection of 
 * candidate teeth for primal separation of simple domino parity inequalities
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_TOOTH_HPP
#define PSEP_TOOTH_HPP

#include "Graph.hpp"

#include <memory>
#include <vector>
#include <map>
#include <unordered_map>

/*
 * If defined, lists of simple teeth will be implemented as unique pointers,
 * else they will be raw pointers.
 */
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
  SimpleTooth() = default;
  
  SimpleTooth(int _root, int _body_start, int _body_end) :
    root(_root), body_start(_body_start), body_end(_body_end),
    slack(INFINITY) {}

  SimpleTooth(int _root, int _body_start, int _body_end, double _slack) :
    root(_root), body_start(_body_start), body_end(_body_end),
    slack(_slack) {}


#ifdef PSEP_TOOTH_UNIQ
  typedef std::unique_ptr<SimpleTooth> Ptr;
#else
  typedef SimpleTooth* Ptr;
#endif

  int root;
  int body_start;
  int body_end;

  double slack;

  int cutgraph_index;

  //returns true iff root lies in the middle of the segment from body_start
  // to body_end
  bool sandwich() const;

  //returns true iff the body segment of the tooth contains node_index
  bool body_contains(const int node_index) const;

  //returns true iff R has the same root as this tooth, and this tooth's
  //body is a subset of R's body
  bool is_subset_of(const SimpleTooth &R) const;
};

/*
 * Structure for representing bodies of simple teeth.
 * start/end are specified relative to a vector of tour nodes, so the segment
 * is
 * tour_nodes[cut_start], tour_nodes[cut_start + 1], ... , tour_nodes[cut_end]
 * slack is the slack of the SEC associated to the segment
 */
struct tooth_seg {
  tooth_seg(int _start, int _end, double _slack) :
    start(_start), end(_end), slack(_slack) {}

  int start;
  int end;
  double slack;
};

/*
 * This class generates and manages a collection of candidate teeth to be 
 * included in a simple domino parity inequality. It is also responsible for
 * nontrivial operations on teeth which require knowledge of the best tour
 * and perm and underlying graph used to define them. 
 */

class CandidateTeeth {
public:
  CandidateTeeth() = default;
  CandidateTeeth(std::vector<int> &_delta, std::vector<int> &_edge_marks,
		 std::vector<int> &_best_tour_nodes,
		 std::vector<int> &_perm,
		 SupportGraph &_G_s,
		 std::vector<int> &_support_indices,
		 std::vector<int> &_support_elist,
		 std::vector<double> &_support_ecap);

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
   * Build a collection of light teeth: teeth for which the slack of the 
   * associated simple tooth inequality is less than 1/2
   */
  int get_light_teeth();

  void print_tooth(const SimpleTooth &T);
  void print_collection();
  
  std::vector<std::vector<SimpleTooth::Ptr>> light_teeth;

  
private:
  void clear_collection();

  /*
   * Reduces the size of the collection of light teeth using the elimination
   * criterion of Lemma 5.5 in Fleischer et al (2006). This reduces the 
   * collection to one of size at most 4m - 2n, without any uncrossing 
   * arguments
   */
  void weak_elim();

  static int add_tooth(std::vector<std::vector<SimpleTooth::Ptr>> &teeth,
		       const int root, const int body_start,
		       const int body_end, const double slack);

  /*
   * this is a callback function to CCcut_linsub_allcuts, u_data should
   * be the LinsubCBData structure below
   */
  static int get_teeth(double cut_val, int cut_start, int cut_end,
		       void *u_data);
  
  std::vector<int> &edge_marks;
  std::vector<int> &best_tour_nodes;
  
  SupportGraph &G_s;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  /*
   * This is the callback handle for CCcut_linsub_allcuts, a pointer to 
   * this should be void-cast and passed to the function along w the callback.
   */
  struct LinsubCBData {
    LinsubCBData(std::vector<std::vector<SimpleTooth::Ptr>> &_cb_teeth,
		 std::vector<int> &_cb_edge_marks,
		 std::vector<int> &_cb_tour_nodes,
		 std::vector<int> &_cb_perm,
		 SupportGraph &_cb_G_s) :
      cb_teeth(_cb_teeth),
      cb_edge_marks(_cb_edge_marks),
      cb_tour_nodes(_cb_tour_nodes), cb_perm(_cb_perm),
      cb_G_s(_cb_G_s) {}


    std::vector<std::vector<SimpleTooth::Ptr>> &cb_teeth;

    std::vector<int> &cb_edge_marks;
    
    std::vector<int> &cb_tour_nodes;
    std::vector<int> &cb_perm;

    SupportGraph &cb_G_s;

    PSEP::tooth_seg *old_seg;

    std::unordered_map<int, double> root_bod_sums;
    std::vector<int> unsorted_roots;
  };

  LinsubCBData cb_data;
  std::vector<int> endmark;
};

}


#endif