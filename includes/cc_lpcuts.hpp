/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief WRAPPERS FOR CONCORDE LP STRUCTURES AND CUT SEPARATORS
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_CC_LPCUTS_HPP
#define PSEP_CC_LPCUTS_HPP

#include "Graph.hpp"
#include "datagroups.hpp"

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
 * @warning CCtsp_lpcut_in is a doubly linked list node, but if you pass one to
 * a CCtsp separation routine, it will store the cuts in a SINGLY LINKED queue.
 */
class LPcutIn {
public:
  LPcutIn();
  ~LPcutIn();

  /** Deletes all non-primal cuts.
   * @pre This function assumes that Concorde heuristic routines were called
   * with the nodes permuted according to the current best tour.
   */
  void filter_primal(PSEP::TourGraph &TG);

  int cut_count() const { return cutcount; }
  bool empty() const { return cutcount == 0; }

  CCtsp_lpcut_in* begin() { return head_cut; } /**< Start of list iterator. */
  CCtsp_lpcut_in* end() { return nullptr; } /**< nullptr for end of list. */
  
  /** Passes address of member pointer for use by separation routines. */
  CCtsp_lpcut_in** pass_ptr() { return &head_cut; }

  /** Passes address of cut count for use by separation routines. */
  int* count_ptr() { return &cutcount; }
  
private:
  void del_cut(CCtsp_lpcut_in *cut);
  
  CCtsp_lpcut_in *head_cut; /**< The raw pointer to the Concorde struct. */
  int cutcount; /**< Number of cuts in the linked list starting at head_cut. */
};

/** Abstract base class for calling Concorde separation routines. 
 * Separator classes based on separation routines from Concorde should derive
 * from this class and provide an implementation of find_cuts which calls
 * an appropriate separation routine or sequence of routines. See below
 * for examples. 
 */
class ConcordeSeparator {
public:
  ConcordeSeparator(PSEP::Data::GraphGroup &_graph_dat,
		    PSEP::Data::BestGroup &_best_dat,
		    PSEP::Data::SupportGroup &_supp_dat,
		    PSEP::TourGraph &_TG, PSEP::Cut::LPcutIn &_cutq) :
    graph_dat(_graph_dat), best_dat(_best_dat), supp_dat(_supp_dat),
    TG(_TG), cutq(_cutq) {}

  /** The call to the separation routine.
   * @returns `true` if cuts are found, `false` otherwise. 
   * @throws std::runtime_error if a call to a concorde routine fails.
   */
  virtual bool find_cuts() = 0;
  
protected:
  
  PSEP::Data::GraphGroup &graph_dat;
  PSEP::Data::BestGroup &best_dat;
  PSEP::Data::SupportGroup &supp_dat;
  
  PSEP::TourGraph &TG;
  
  PSEP::Cut::LPcutIn &cutq;
};

/** Exact separation of segment cut subtours. */
class SegmentCuts : public ConcordeSeparator {
public:
  SegmentCuts(PSEP::Data::GraphGroup &g_dat, PSEP::Data::BestGroup &b_dat,
	      PSEP::Data::SupportGroup &s_dat, PSEP::TourGraph &TG,
	      PSEP::Cut::LPcutIn &cutq) :
    ConcordeSeparator(g_dat, b_dat, s_dat, TG, cutq) {}

  bool find_cuts();
};

/** Primal separation of comb ineqalities via standard block comb heuristic. */
class BlockCombs : public ConcordeSeparator {
public:
  BlockCombs(PSEP::Data::GraphGroup &g_dat, PSEP::Data::BestGroup &b_dat,
	     PSEP::Data::SupportGroup &s_dat, PSEP::TourGraph &TG,
	     PSEP::Cut::LPcutIn &cutq) :
    ConcordeSeparator(g_dat, b_dat, s_dat, TG, cutq) {}

  bool find_cuts();
};

/** Primal separation of blossoms via standard fast blossom heuristics. 
 * Calls Padberg-Hong odd component blossoms, and then Grotschel-Holland fast
 * blossoms if none are found.
 */
class FastBlossoms : public ConcordeSeparator {
public:
  FastBlossoms(PSEP::Data::GraphGroup &g_dat, PSEP::Data::BestGroup &b_dat,
	       PSEP::Data::SupportGroup &s_dat, PSEP::TourGraph &TG,
	       PSEP::Cut::LPcutIn &cutq) :
    ConcordeSeparator(g_dat, b_dat, s_dat, TG, cutq) {}

  bool find_cuts();
};

}
}

#endif
