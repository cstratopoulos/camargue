/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Wrappers for Concorde LP structures and separators.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CC_LPCUTS_HPP
#define CMR_CC_LPCUTS_HPP

#include "graph.hpp"
#include "datagroups.hpp"

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

#include <memory>
#include <vector>


namespace CMR {

/** Classes and functions related to cut separation. */
namespace Sep {

/** Wrapper to the Concorde CCtsp_lpgraph structure.
 * This class constructs a CCtsp_lpgraph which corresponds to a tour specified
 * by the constructor arguments. It is used to check whether cuts found by
 * Concorde standard heuristics are tight at the current tour.
 */
class TourGraph {
public:
    TourGraph() noexcept;

    TourGraph(const std::vector<double> &tour_edges,
              const std::vector<Graph::Edge> &edges,
              const std::vector<int> &perm);

    TourGraph(TourGraph &&T) noexcept;
    TourGraph &operator=(TourGraph &&T) noexcept;

    ~TourGraph();

    TourGraph(const TourGraph &T) = delete;
    TourGraph &operator=(const TourGraph &T) = delete;



    CCtsp_lpgraph* pass_ptr() { return &L; }
    double* tour_array() { return &d_tour[0]; }
    int node_count() const { return L.ncount; }

private:
  std::vector<double> d_tour;
  CCtsp_lpgraph L;
};

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

    CCtsp_lpcut_in *begin() { return head_cut.get(); }
    const CCtsp_lpcut_in *begin() const { return head_cut.get(); }

    void filter_primal(TourGraph &TG);

    void clear();

private:
  struct hungry_delete {
    void operator()(CCtsp_lpcut_in *cut) const {
      for(auto it = cut; it; it = cut){
	cut = it->next;
	CCtsp_free_lpcut_in(it);
	CC_IFFREE(it, CCtsp_lpcut_in);
      }
    }
  };
  std::unique_ptr<CCtsp_lpcut_in, hungry_delete> head_cut;
  int cutcount;
};


/** Abstract base class for calling Concorde separation routines.
 * Separator classes based on separation routines from Concorde should derive
 * from this class and provide an implementation of find_cuts which calls
 * an appropriate separation routine or sequence of routines. See FastBlossoms,
 * BlockCombs, or SegmentCuts for examples.
 */
class ConcordeSeparator {
public:
    ConcordeSeparator(std::vector<int> &supp_elist,
                      std::vector<double> &supp_ecap,
                      TourGraph &_TG, Sep::LPcutList &_cutq) :
        elist(supp_elist), ecap(supp_ecap), TG(_TG), cutq(_cutq) {}

    /** The call to the separation routine.
     * @returns `true` if cuts are found, `false` otherwise.
     * @throws std::runtime_error if a call to a concorde routine fails.
     */
    virtual bool find_cuts() = 0;

protected:
    std::vector<int> &elist;
    std::vector<double> &ecap;

    TourGraph &TG;

    Sep::LPcutList &cutq;
};

/** Exact separation of segment cut subtours. */
class SegmentCuts : ConcordeSeparator {
public:
    SegmentCuts(std::vector<int> &elist, std::vector<double> &ecap,
                TourGraph &TG, Sep::LPcutList &cutq) :
        ConcordeSeparator(elist, ecap, TG, cutq) {}

    /** Finds subtours arising from intervals of the current best tour. */
    bool find_cuts();
};

/** Standard separation of connect cuts. */
class ConnectCuts : ConcordeSeparator {
public:
    ConnectCuts(std::vector<int> &elist, std::vector<double> &ecap,
                TourGraph &TG, Sep::LPcutList &cutq) :
        ConcordeSeparator(elist, ecap, TG, cutq) {}

    /** Finds subtours arising from connected components. */
    bool find_cuts();
};

/** Primal separation of comb ineqalities via standard block comb heuristic. */
class BlockCombs : ConcordeSeparator {
public:
    BlockCombs(std::vector<int> &elist, std::vector<double> &ecap,
               TourGraph &TG, Sep::LPcutList &cutq) :
        ConcordeSeparator(elist, ecap, TG, cutq) {}

    /** Returns true if block combs are found and some are tight at tour. */
    bool find_cuts();
};

/** Primal separation of blossoms via standard fast blossom heuristics. */
class FastBlossoms : ConcordeSeparator {
public:
    FastBlossoms(std::vector<int> &elist, std::vector<double> &ecap,
                 TourGraph &TG, Sep::LPcutList &cutq) :
        ConcordeSeparator(elist, ecap, TG, cutq) {}

    /** Returns true if blossoms are found and some are tight at best tour.
     * First calls the Padberg-Hong odd component blossom heuristic, and then
     * the Grotschell-Holland heuristic if no odd component blossoms are found.
     */
    bool find_cuts();
};

}
}

#endif
