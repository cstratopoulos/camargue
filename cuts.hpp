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
#include "tooth.hpp"
#include "datagroups.hpp"
#include "Graph.hpp"

namespace PSEP {
constexpr int max_add = 4;
constexpr int seg_q_max = max_add;
constexpr int blossom_q_max = 15;

/*
 * Structure for storing segment cuts: subtour inequalities arising from
 * segments of the current best tour.
 * start, end - the start and endpoints of the segment, indicated as 
 *              indices for accessing the vector best_tour_nodes. 
 *              Thus the associated segment is best_tour_nodes[start] to
 *              best_tour_nodes[end]
 * viol - the amount by which the cut is violated
 */
struct seg {
  seg(int _start, int _end, double _viol) :
    start(_start), end(_end), viol(_viol) {}

  int start;
  int end;
  double viol;
};

/*
 * Structure for storing blossom inequalities, aka 2-matching inequalities.
 * These are simple comb inequalities where each tooth is an edge.
 */
struct blossom {
  blossom(std::vector<int> &_handle, int _cut_edge, double _val) :
    handle(_handle), cut_edge(_cut_edge), cut_val(_val){}

  std::vector<int> handle;

  /* cut_edge represents the edge {u, v} for which a minimum uv-cut is 
   * computed in the separation routine. From e, handle, and the 
   * vector best_tour_edges it is possible to extract the blossom inequality
   * See Letchford-Lodi Primal separation algorithms for details
   */
  int cut_edge;
  double cut_val;
};

/* The cuts below are dummy structures which are not actually used at the
 * moment but may be useful if cut pools or column gen are implemented */
  
/* Struct for storing simple DP inequalities */
struct domino {
  domino(){}
  domino(std::vector<int> &_handle,
	 std::vector<std::shared_ptr<CandTooth::SimpleTooth>> _teeth) :
    handle(_handle), used_teeth(_teeth) {}

  std::vector<int> handle;
  std::vector<std::shared_ptr<CandTooth::SimpleTooth>> used_teeth;
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
  virtual int cut_call() = 0;

protected:
  virtual int separate() = 0;
  virtual int add_cut() = 0;
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
  CutQueue(const int cap) : q_capacity(cap) {}
  
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
  
  void pop_front() { cut_q.pop_front(); }  
  bool empty() const { return cut_q.empty(); }
  int size() const { return cut_q.size(); }

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
			std::vector<int> &rmatind, std::vector<double> &rmatval,
			char &sense, double &rhs);

private:
  std::vector<Edge> &edges;
  std::vector<int> &delta;
  std::vector<int> &edge_marks;
  IntPairMap &edge_lookup;
};

}

/*
 *       FORWARD DECLARATIONS OF PARTIAL TEMPLATE SPECIALIZATIONS
 */

template<>
void PSEP::CutQueue<PSEP::HyperGraph>::push_front(const PSEP::HyperGraph &H);

template<>
void PSEP::CutQueue<PSEP::HyperGraph>::pop_front();

#endif
