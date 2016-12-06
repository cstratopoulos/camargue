/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief WRAPPERS FOR CONCORDE LP STRUCTURES AND CUT SEPARATORS
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CC_LPCUTS_HPP
#define CMR_CC_LPCUTS_HPP

#include "Graph.hpp"
#include "datagroups.hpp"

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

#include <memory>

namespace std {

template<>
struct default_delete<CCtsp_lpcut_in> {
  void operator()(CCtsp_lpcut_in *cut) const {
    for(auto it = cut; it; it = cut){
      cut = it->next;
      CCtsp_free_lpcut_in(it);
      CC_IFFREE(it, CCtsp_lpcut_in);
    }
  }
};

}


namespace CMR {

/** Classes and functions related to cut separation. */
namespace Cut {

class LPcutList {
public:
  LPcutList() noexcept;
  LPcutList(CCtsp_lpcut_in *head, int count) noexcept;
  LPcutList(LPcutList &&L) noexcept;
  LPcutList(const LPcutList &L) = delete;

  LPcutList &operator=(LPcutList &&L) noexcept;
  LPcutList &operator=(const LPcutList &L) = delete;

  int size() const { return cutcount; }
  bool empty() const { return cutcount == 0; }
  
  CCtsp_lpcut_in* begin() { return head_cut.get(); }

  void filter_primal(CMR::TourGraph &TG);
  
private:  
  std::unique_ptr<CCtsp_lpcut_in> head_cut;
  int cutcount;
};


/** Wrapper for Concorde CCtsp_lpcut_in structure. 
 * This class is meant to provide an extremely limited (but memory safe!) 
 * interface, for no other purpose but to call certain cut separation routines
 * provided by Concorde. It models unique ownership over a list of cuts. Copy
 * operations are deleted, and move operations transfer ownership.
 * @warning CCtsp_lpcut_in is a doubly linked list node, but if you pass one to
 * a CCtsp separation routine, it will store the cuts in a SINGLY LINKED queue.
 */
class LPcutIn {
public:
  LPcutIn() noexcept; /**< Default construct an empty list with null head. */
  
  LPcutIn(const LPcutIn &C) = delete; /**< Deleted copy constructor. */
  
  LPcutIn(LPcutIn &&C) noexcept; /**< Construct, transferring ownership. */

  LPcutIn &operator=(const LPcutIn &C) = delete; /**< Deleted copy assign. */

  LPcutIn &operator=(LPcutIn &&C) noexcept; /**< Assign, transfer ownership.*/
  
  ~LPcutIn(); /**< Traverse the list, freeing every entry. */

  /** Deletes all non-primal cuts.
   * @pre This function assumes that Concorde heuristic routines were called
   * with the nodes permuted according to the current best tour.
   */
  void filter_primal(CMR::TourGraph &TG);

  int cut_count() const { return cutcount; }
  bool empty() const { return cutcount == 0; }

  CCtsp_lpcut_in* begin() { return head_cut; } /**< Start of list iterator. */
  CCtsp_lpcut_in* end() { return (CCtsp_lpcut_in *) NULL; } /**< List end. */
  
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
 * an appropriate separation routine or sequence of routines. See FastBlossoms,
 * BlockCombs, or SegmentCuts for examples. 
 */
class ConcordeSeparator {
public:
  ConcordeSeparator(CMR::Data::GraphGroup &_graph_dat,
		    CMR::Data::BestGroup &_best_dat,
		    CMR::Data::SupportGroup &_supp_dat,
		    CMR::TourGraph &_TG, CMR::Cut::LPcutList &_cutq) :
    graph_dat(_graph_dat), best_dat(_best_dat), supp_dat(_supp_dat),
    TG(_TG), cutq(_cutq) {}

  /** The call to the separation routine.
   * @returns `true` if cuts are found, `false` otherwise. 
   * @throws std::runtime_error if a call to a concorde routine fails.
   */
  virtual bool find_cuts() = 0;
  
protected:
  
  CMR::Data::GraphGroup &graph_dat;
  CMR::Data::BestGroup &best_dat;
  CMR::Data::SupportGroup &supp_dat;
  
  CMR::TourGraph &TG;
  
  CMR::Cut::LPcutList &cutq;
};

/** Exact separation of segment cut subtours. */
class SegmentCuts : public ConcordeSeparator {
public:
  SegmentCuts(CMR::Data::GraphGroup &g_dat, CMR::Data::BestGroup &b_dat,
	      CMR::Data::SupportGroup &s_dat, CMR::TourGraph &TG,
	      CMR::Cut::LPcutList &cutq) :
    ConcordeSeparator(g_dat, b_dat, s_dat, TG, cutq) {}

  /** Finds subtours arising from intervals of the current best tour. */
  bool find_cuts(); 
};

/** Primal separation of comb ineqalities via standard block comb heuristic. */
class BlockCombs : public ConcordeSeparator {
public:
  BlockCombs(CMR::Data::GraphGroup &g_dat, CMR::Data::BestGroup &b_dat,
	     CMR::Data::SupportGroup &s_dat, CMR::TourGraph &TG,
	     CMR::Cut::LPcutList &cutq) :
    ConcordeSeparator(g_dat, b_dat, s_dat, TG, cutq) {}

  /** Returns true if block combs are found and some are tight at best tour. */
  bool find_cuts();
};

/** Primal separation of blossoms via standard fast blossom heuristics. */
class FastBlossoms : public ConcordeSeparator {
public:
  FastBlossoms(CMR::Data::GraphGroup &g_dat, CMR::Data::BestGroup &b_dat,
	       CMR::Data::SupportGroup &s_dat, CMR::TourGraph &TG,
	       CMR::Cut::LPcutList &cutq) :
    ConcordeSeparator(g_dat, b_dat, s_dat, TG, cutq) {}

  /** Returns true if blossoms are found and some are tight at best tour.
   * First calls the Padberg-Hong odd component blossom heuristic, and then
   * the Grotschell-Holland heuristic if no odd component blossoms are found.
   */
  bool find_cuts();
};

}
}

#endif
