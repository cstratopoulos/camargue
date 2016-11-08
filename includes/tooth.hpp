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

/** Structure for storing candidate simple teeth.
 * This structure holds a candidate simple tooth; a collection of these may
 * potentially comprise a simple domino parity inequalities. A simple tooth is
 * a root node, \p root, and a segment of a
 * best known tour, stored as best tour node indices. Hence the root 
 * node is `best_tour_nodes[root]`, and the segment is
 * `best_tour_nodes[body_start], best_tour_nodes[body_start + 1], ...
 *                                best_tour_nodes[body_end]`,
 * incrementing indices mod the number of nodes in the graph.
 * A SimpleTooth is meaningless without management by its CandidateTooth
 * class, which is responsible for all nontrivial operations on the tooth. 
 */

struct SimpleTooth {
  SimpleTooth() = default;

  /** Constructs a SimpleTooth with root and body segment. */
  SimpleTooth(int _root, int _body_start, int _body_end) :
    root(_root), body_start(_body_start), body_end(_body_end),
    slack(INFINITY) {}

  /** Constructs a SimpleTooth with root, body segment, and slack. */
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

  /** The slack of the simple tooth inequality. 
   * Suppose the tooth has \p root \f$ i \f$, and body \f$ S \f$ equal to 
   * `best_tour_nodes[body_start]`, `best_tour_nodes[body_start + 1]`, ...,
   * `best_tour_nodes[body_end]`.
   * If \f$ x^* \f$ is the %LP solution, then \p slack is \f[
   * 2x(E(S)) + x(E(i:S)) \le 2 |S| - 1. \f]
   */
  double slack;

  int cutgraph_index; /** Labelling index for use in DPCutGraph. */

  //returns true iff root lies in the middle of the segment from body_start
  // to body_end
  bool sandwich() const;

  /** Returns true iff the body segment of the tooth contains \p node_index. */
  bool body_contains(const int node_index) const;

  /** Returns true if this tooth's body is a subset of \p R's body. */
  bool is_subset_of(const SimpleTooth &R) const;
};

/** Structure for representing bodies of simple teeth.
 * \p start and \p end are specified relative to a vector of tour nodes, 
 * so the segment is
 * `tour_nodes[cut_start], tour_nodes[cut_start + 1], ... , tour_nodes[cut_end]`
 */
struct tooth_seg {
  /** Construct a segment with given \p start, \p end, and \p slack. */
  tooth_seg(int _start, int _end, double _slack) :
    start(_start), end(_end), slack(_slack) {}

  int start;
  int end;
  double slack; /** The slack of the SEC associated to the segment. */
};

/** Class for generating and managing a collection of candidate SimpleTeeth.
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

  /** The number of nodes in \p T's body. */
  int body_size(const SimpleTooth &T);

  /** The relative complement of a SimpleTooth.
   * If \p T has root \f$ i \f$ and body \f$ S \f$, then `complement(T)` has
   * root \f$ i \f$ and body \f$ V \setminus (S cup \left\{ i \right\}) \f$.
   */
  void complement(SimpleTooth &T);

  /** Is the body of \p T a subset of the body of \p R.
   * For teeth \p T, \p R with the same \p root, returns whether the body of
   * \p T is a subset of the body of \p R, storing the result in \p result.
   * @returns 0 if success, 1 if failure.
   * @todo generalize to arbitrary teeth.
   */
  int body_subset(const SimpleTooth &T, const SimpleTooth &R, bool &result);

  /** Build a collection of light simple teeth. */
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
