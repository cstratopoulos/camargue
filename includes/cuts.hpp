/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                            
 *                CUT TEMPLATE AND CUT STRUCTURE DEFINITIONS          
 *
 * This file contains some structure definitions for types of cuts used in 
 * the TSP solver, as well as a Cut (separator) template and a CutQueue template
 * for use by the Cutcontrol class
 *                                                                            
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PSEP_CUTS_H
#define PSEP_CUTS_H

#include <memory>
#include <vector>

#include "PSEP_util.hpp"
#include "setbank.hpp"
#include "datagroups.hpp"
#include "Graph.hpp"
#include "tooth.hpp"

namespace PSEP {

/** Structure for storing segment cuts.
 * Segment cuts are subtour inequalities arising from segments of the current
 *  best tour. If `tour` is a vector of tour nodes, the segment is given by
 * `tour[start], tour[start + 1], ..., tour[end]`. If \f$ x^* \f$ is an lp 
 * solution and \f$ S \f$ is the segment, \p cutval is \f$ x^*(\delta(S)) \f$ .
 */
struct seg {
  seg() = default;
  seg(int _start, int _end, double _cutval) :
    start(_start), end(_end), cutval(_cutval) {}

  int start;
  int end;
  double cutval;
};

/** Structure for storing blossom inequalities from exact primal separation.
 * This structure stores blossom inequalities found by the exact primal
 * separation algorithm of Letchford-Lodi in Primal Separation Algorithms
 * Revisited. 
 */
struct blossom {
  blossom(std::vector<int> &_handle, int _cut_edge, double _val) :
    handle(_handle), cut_edge(_cut_edge), cut_val(_val){}

  
  std::vector<int> handle;
  
  int cut_edge;

  double cut_val;
};

/** Structure for storing blossoms found by fast heuristics. */
struct fastblossom {
  fastblossom(std::vector<int> &_handle, std::vector<int> &_edge_indices) :
    handle(_handle), edge_indices(_edge_indices) {}

  std::vector<int> handle;
  std::vector<int> edge_indices;
};

/** Structure for storing simple DP inequalities.
 * In all cases, all numbers and indices refer to position in some `tour_nodes`
 * vector. They must be translated by deferencing, e.g., if `(i, j)` is an
 * edge in \p nonneg_edges, it refers to the edge 
 * `tour_nodes[i], tour_nodes[j]`.
 */
struct dominoparity {
  dominoparity() = default;
  dominoparity(std::vector<PSEP::SimpleTooth*> &_used_teeth,
	       std::vector<int> &_degree_nodes,
	       std::vector<IntPair> &_nonneg_edges) :
    used_teeth(_used_teeth), degree_nodes(_degree_nodes),
    nonneg_edges(_nonneg_edges) {}

  /** Simple tooth inequalities to be aggregated to get the simple DP ineq. */
  std::vector<PSEP::SimpleTooth*> used_teeth;

  /** The handle of the inequality, nodes where the degree eqns are used. */
  std::vector<int> degree_nodes;

  /** Edges \f$ e \f$ for which \f$ x_e \ge 0 \f$ is used. */
  std::vector<IntPair> nonneg_edges;
};

/* cut mimicking the parameters used to add a row in CPLEX; used for safe
 * Gomory cut separation
 */
struct safeGMI {
  safeGMI(std::vector<int> &_rmatind, std::vector<double> &_rmatval,
	  char _sense, double _RHS) :
    rmatind(_rmatind), rmatval(_rmatval), sense(_sense), RHS(_RHS) {}

  std::vector<int> rmatind;
  std::vector<double> rmatval;
  char sense;
  double RHS;
};


/*
 * This pure abstract class defines the interface to a separation routine
 * for a given cut type
 *
 * cut_call - the wrapper function calling all the protected methods
 * separate - invokes the separation routine for cuts of type cut_t
 * add_cut - the function for actually adding the row
 *           TODO: will just add the hypergraph to queue??
 * see segments.hpp, blossoms.hpp, dominos.hpp, etc for examples
 */
template<typename cut_t> class Cut {
public:
  virtual int cutcall();// = 0;

protected:
  virtual int separate();// = 0;
  virtual int add_cuts();// = 0;
};


/*
 * This template class provides an interface for dealing with a queue of
 * representations of cuts, cut_rep. In practice this 
 * will be a HyperGraph (see setbank.hpp), or a structure for storing simple
 * DP inequalities (see dominos.hpp, simpleDP.hpp), or the locally used 
 * structure from within a separator (see segments.cpp, blossoms.cpp)
 */
  
template<typename cut_rep>
class CutQueue {
public:
  CutQueue(const int cap) : q_capacity(cap), q_fresh(true) {}
  
  //the max number of cuts to be stored in the queue
  const int q_capacity;

  //returns a reference to the most recently added cut
  const cut_rep &peek_front() const { return cut_q.front(); }
  
  //pushes a new cut to the front of the queue, popping the oldest one from
  //the back if we are at capacity
  void push_front(const cut_rep &H)
  {
    cut_q.push_front(H);
    if(cut_q.size() > q_capacity) cut_q.pop_back();
  }

  void push_back(const cut_rep &H)
  {
    if(cut_q.size() >= q_capacity) cut_q.pop_back();
    cut_q.push_back(H);
  }
  
  void pop_front() { cut_q.pop_front(); }  
  bool empty() const { return cut_q.empty(); }
  int size() const { return cut_q.size(); }

  bool q_fresh;

private:
  std::list<cut_rep> cut_q;
};

class CutTranslate {
public:
  CutTranslate(Data::GraphGroup &GraphGroup) :
    edges(GraphGroup.m_graph.edges),
    delta(GraphGroup.delta),
    edge_marks(GraphGroup.edge_marks),
    edge_lookup(GraphGroup.m_graph.edge_lookup) {}

  int get_sparse_row(const HyperGraph &H, std::vector<int> &rmatind,
		     std::vector<double> &rmatval, char &sense, double &rhs);
  int get_sparse_row_if(bool &violated, const HyperGraph &H,
			const std::vector<double> &x,
			std::vector<int> &rmatind, std::vector<double> &rmatval,
			char &sense, double &rhs);
  int is_cut_violated(bool &violated, const HyperGraph &H,
		      std::vector<double> &x);

private:
  template<typename number_t>
  void get_activity(double &activity, const std::vector<number_t> &x,
		    const std::vector<int> &rmatind,
		    const std::vector<double> &rmatval)
  {
    activity = 0;
    for(int i = 0; i < rmatind.size(); i++){
      int index = rmatind[i];
      activity += x[index] * rmatval[i];
    }
  }
  
  std::vector<Edge> &edges;
  std::vector<int> &delta;
  std::vector<int> &edge_marks;
  IntPairMap &edge_lookup;
};

/*
 *       FORWARD DECLARATIONS OF PARTIAL TEMPLATE SPECIALIZATIONS
 */


template<>
void CutQueue<HyperGraph>::push_front(const HyperGraph &H);

template<>
void CutQueue<HyperGraph>::push_back(const HyperGraph &H);

template<>
void CutQueue<HyperGraph>::pop_front();
  
}



#endif
