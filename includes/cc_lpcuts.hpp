/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Wrappers for Concorde cut structures/separators.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CC_LPCUTS_H
#define CMR_CC_LPCUTS_H

#include "graph.hpp"
#include "datagroups.hpp"

#include <stdexcept>
#include <string>
#include <memory>
#include <vector>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}


namespace CMR {

/// Classes and functions related to cut separation.
namespace Sep {

/// Wrapper to Concorde CCtsp_lpgraph for pricing cuts at tours.
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

/// Management of Concorde lpcut_in linked list.
class LPcutList {
public:
    LPcutList() noexcept;

    /// Construct a list of \p count cuts with head cut \p head.
    LPcutList(CCtsp_lpcut_in *head, int count) noexcept;

    LPcutList(LPcutList &&L) noexcept;
    LPcutList &operator=(LPcutList &&L) noexcept;

    void push_front(CCtsp_lpcut_in *new_head);

    void splice(LPcutList &&L);

    int size() const { return cutcount; }
    bool empty() const { return cutcount == 0; }

    CCtsp_lpcut_in *begin() { return head_cut.get(); }
    const CCtsp_lpcut_in *begin() const { return head_cut.get(); }

    void filter_primal(TourGraph &TG);

    void clear();

private:
    /// Deleter to clear the linked list.
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


/// Abstract base class for calling Concorde separation routines.
/// See ConcordeSeparator and type aliases for sample usage, or LocalCuts.
class CCsepBase {
public:
    CCsepBase(std::vector<int> &supp_elist,
              std::vector<double> &supp_ecap,
              TourGraph &_TG, Sep::LPcutList &_cutq) :
        elist(supp_elist), ecap(supp_ecap), TG(_TG), cutq(_cutq) {}

    /// Call the separation routine, returning true iff cuts are found.
    virtual bool find_cuts() = 0;

    bool filter_primal = true; //!< Should only tight cuts be kept.

protected:
    std::vector<int> &elist;
    std::vector<double> &ecap;

    TourGraph &TG;

    Sep::LPcutList &cutq;
};

/// Primal separation of non-template local cuts via standard heuristics.
class LocalCuts : public CCsepBase {
public:
    LocalCuts(std::vector<int> &elist, std::vector<double> &ecap,
              TourGraph &TG, Sep::LPcutList &cutq, int seed) :
        CCsepBase(elist, ecap, TG, cutq), random_seed(seed) {}

    bool find_cuts(); //!< Returns true if tight local cuts are found.

    static constexpr int MaxChunkSize = 16;
    int current_max = 8;
    int random_seed;
    bool spheres = false;
};

/// Alias declaration for Concorde separator call prototype.
using CCsepCall = decltype(&CCtsp_segment_cuts);

/** Class template for straightforward Concorde separation routines.
 * @tparam sep_fn the Concorde function to call.
 * @tparam check_filter_primal if false, the value CCsepBase#filter_primal
 * will be disgregarded: no primal filtering ever takes place. Useful for
 * SegmentCuts, where the separation routine is already primal, or for other
 * SEC routines where we only invoke them if we want non-tight cuts.
 * @tparam fn_name the name of the function, in case an error message needs
 * to be printed.
 * This class template can be used to add Concorde separators to Camargue
 * in a simple way. See type aliases SegmentCuts, FastBlossoms, etc. for sample
 * usage.
 */
template <CCsepCall sep_fn, bool check_filter_primal, const char *fn_name>
class ConcordeSeparator : public CCsepBase {
public:
    ConcordeSeparator(std::vector<int> &elist, std::vector<double> &ecap,
                      TourGraph &TG, Sep::LPcutList &cutq) :
        CCsepBase(elist, ecap, TG, cutq) {}

    bool find_cuts()
        {
            int cutcount = 0;
            CCtsp_lpcut_in *head = NULL;
            std::string err_msg(fn_name); err_msg += " failed";


            if (sep_fn(&head, &cutcount, TG.node_count(), ecap.size(),
                       &elist[0], &ecap[0]))
                throw std::runtime_error(err_msg);

            if (cutcount == 0)
                return false;

            cutq = LPcutList(head, cutcount);

            if (check_filter_primal && filter_primal)
                cutq.filter_primal(TG);

            return !cutq.empty();
        }
};

/**@name String names for ConcordeSeparator template instantiation.
 * These need to be defined to pass a string literal as a template parameter,
 * see eg http://www.comeaucomputing.com/techtalk/templates/#stringliteral
 */
///@{
constexpr char seg_fname[] = "CCtsp_segment_cuts";
constexpr char con_fname[] = "CCtsp_connect_cuts";
constexpr char ex_fname[] = "CCtsp_exact_subtours";
constexpr char f2m_fname[] = "CCtsp_fastblossom";
constexpr char gh2m_fname[] = "CCtsp_ghfastblossom";
constexpr char blk_fname[] = "CCtsp_block_combs";
///@}

/**@name ConcordeSeparator explicit instantiations.
 * Except for SegmentCuts, all of these are standard heuristics.
 */
///@{

/// Exact primal SEC separation.
using SegmentCuts = ConcordeSeparator<CCtsp_segment_cuts, false, seg_fname>;

/// Standard connected component SEC generation.
using ConnectCuts = ConcordeSeparator<CCtsp_connect_cuts, false, con_fname>;

/// Exact standard SEC separation.
using ExactSub = ConcordeSeparator<CCtsp_exact_subtours, false, ex_fname>;

/// Odd component fast blossoms.
using FastBlossoms = ConcordeSeparator<CCtsp_fastblossom, true, f2m_fname>;

/// Gr\"otschel-Holland blossom heuristic.
using GHblossoms = ConcordeSeparator<CCtsp_ghfastblossom, true, gh2m_fname>;

/// Wrapper function because CCtsp_block_combs takes a silent parameter.
/// @remark All (?) sep routines have this prototype in later versions of
/// Concorde.
inline int BlkCombCall(CCtsp_lpcut_in **c, int *cutcount,
                       int ncount, int ecount, int *elist, double *ecap)
{
    return CCtsp_block_combs(c, cutcount, ncount, ecount, elist, ecap, 1);
}

/// Block comb separation.
using BlockCombs = ConcordeSeparator<BlkCombCall, true, blk_fname>;

///@}

}
}

#endif
