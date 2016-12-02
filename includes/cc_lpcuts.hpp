/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief WRAPPERS FOR CONCORDE LP STRUCTURES
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_CC_LPCUTS_HPP
#define PSEP_CC_LPCUTS_HPP

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace PSEP {

/** Classes and functions related to cut separation. */
namespace Cut {

/** Wrapper for Concorde CCtsp_lpcut_in structure. 
 * This class is meant to provide an extremely limited (but memory safe!) 
 * interface, for no other purpose but to call certain cut separation routines
 * provided by Concorde. 
 * @
 */
class CCwrapper {
public:
  CCwrapper();
  ~CCwrapper();

  /** Deletes all non-primal cuts.
   * @pre This function assumes that Concorde heuristic routines were called
   * with the nodes permuted according to the current best tour.
   */
  void filter_primal(); 

  int cut_count() const { return cutcount; }
  bool empty() const { return cutcount == 0; }
  
  /** Passes address of member pointer for use by separation routines. */
  CCtsp_lpcut_in** pass_ptr() { return &cc_cut; }

  /** Passes address of cut count for use by separation routines. */
  int* count_ptr() { return &cutcount; }
  
private:
  CCtsp_lpcut_in *cc_cut; /**< The raw pointer to the Concorde struct. */
  int cutcount; /**< Number of cuts in the linked list starting at cc_cut. */
};


}
}

#endif
