/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                            
 *                CUT TEMPLATE AND CUT STRUCTURE DEFINITIONS          
 *
 * This file contains some structure definitions for types of cuts used in 
 * the TSP solver, as well as a dummy cut template for use by the CutControl
 * class
 *                                                                            
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PSEP_CUTS_H
#define PSEP_CUTS_H

#include <memory>
#include <vector>

#include "tooth.h"

namespace PSEP {

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
   * see segments.h, blossoms.h, dominos.h, etc for examples
   */
  template<typename cut_t> class Cut {
  public:
    virtual int cut_call() = 0;

  protected:
    virtual int separate() = 0;
    virtual int add_cut() = 0;
  };
}

#endif
