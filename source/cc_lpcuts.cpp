#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>
#include <utility>

extern "C" {
#include <concorde/INCLUDE/localcut.h>
}

using std::vector;

using std::cout;
using std::endl;
using std::cerr;

using std::runtime_error;

using lpcut_in = CCtsp_lpcut_in;

namespace CMR {
namespace Sep {

TourGraph::TourGraph() noexcept { CCtsp_init_lpgraph_struct(&L); }

TourGraph::TourGraph(const vector<double> &tour_edges,
		     const vector<Graph::Edge> &edges,
                     const vector<int> &perm) try
    : d_tour(tour_edges)
{
    vector<int> elist;
    int ncount = perm.size();
    int ecount = edges.size();

    for (const Graph::Edge &e : edges) {
        elist.push_back(perm[e.end[0]]);
        elist.push_back(perm[e.end[1]]);
    }

  CCtsp_init_lpgraph_struct(&L);

  if (CCtsp_build_lpgraph(&L, ncount, ecount, &elist[0], (int *) NULL))
      throw runtime_error("CCtsp_build_lpgraph failed.");

  if (CCtsp_build_lpadj(&L, 0, ecount))
      throw runtime_error("CCtsp_build_lpadj failed.");

} catch (const std::exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("TourGraph constructor failed.");
}

TourGraph::TourGraph(TourGraph &&T) noexcept : d_tour(std::move(T.d_tour))
{
    CCtsp_free_lpgraph(&L);
    L = T.L;

    CCtsp_init_lpgraph_struct(&T.L);
    T.d_tour.clear();
}

TourGraph &TourGraph::operator=(TourGraph &&T) noexcept
{
    d_tour = std::move(T.d_tour);

    CCtsp_free_lpgraph(&L);
    L = T.L;

    CCtsp_init_lpgraph_struct(&T.L);
    T.d_tour.clear();

    return *this;
}

TourGraph::~TourGraph() {  CCtsp_free_lpgraph(&L); }


LPcutList::LPcutList() noexcept : head_cut(), cutcount(0) {}

LPcutList::LPcutList(CCtsp_lpcut_in *head, int count) noexcept :
  head_cut(head), cutcount(count) {}

LPcutList::LPcutList(LPcutList &&L) noexcept :
  head_cut(std::move(L.head_cut)), cutcount(L.cutcount) { L.cutcount = 0; }

LPcutList &LPcutList::operator=(LPcutList &&L) noexcept {
  head_cut = std::move(L.head_cut);
  cutcount = L.cutcount;
  L.cutcount = 0;
  return *this;
}

void LPcutList::push(CCtsp_lpcut_in *new_head)
{
    ++cutcount;
    new_head->next = head_cut.release();
    head_cut.reset(new_head);
}

void LPcutList::pop()
{
    if (empty() || !head_cut)
        throw runtime_error("Popped from empty LPcutList");

    CCtsp_lpcut_in *new_head = head_cut->next;
    head_cut->next = nullptr;
    head_cut.reset(new_head);

    --cutcount;
}

void LPcutList::splice(LPcutList &&L)
{
    if (empty()) {
        head_cut = std::move(L.head_cut);
        cutcount = L.cutcount;
    } else {
        auto it = begin();
        while(it->next != nullptr)
            it = it->next;
        it->next = L.head_cut.release();
        cutcount += L.cutcount;
    }

    L.cutcount = 0;
}

void LPcutList::filter_primal(TourGraph &TG)
{
    if (cutcount == 0 || !head_cut) return;

    lpcut_in *current = head_cut.get();
    lpcut_in *prev = current;

    while (current) {
        double slack = CCtsp_cutprice(TG.pass_ptr(), current, TG.tour_array());

        if (slack != 0) {
            --cutcount;

            if (current == head_cut.get()) {
                current = head_cut->next;
                head_cut->next = nullptr;
                head_cut.reset(current);
                current = head_cut.get();
                prev = current;
            } else {
                prev->next = current->next;
                CCtsp_free_lpcut_in(current);
                CC_IFFREE(current, lpcut_in);
                current = prev->next;
            }
        } else {
            prev = current;
            current = current->next;
        }
    }
}

void LPcutList::clear()
{
    cutcount = 0;
    head_cut.reset();
}



bool LocalCuts::find_cuts()
{
    CCchunk_flag flags;
    CCchunk_localcut_timer lc_timer;
    CCchunk_init_localcut_timer(&lc_timer);

    flags.dummy = 0;
    flags.permute = 0;
    flags.weighted = 0;
    flags.spheres = spheres;
    flags.uncivilized = 0;
    flags.noshrink = 0;
    flags.nolift = 0;

    int cutcount = 0;
    lpcut_in *head = NULL;

    CCrandstate rstate;
    CCutil_sprand(random_seed, &rstate);

    flags.maxchunksize = current_max;
    flags.spheresize   = current_max - 2;

    if (CCchunk_localcuts(&head, &cutcount, TG.node_count(), ecap.size(),
                          &elist[0], &ecap[0], 0.0, flags, &lc_timer, 1,
                          &rstate))
        throw runtime_error("CCchunk_localcuts failed.");

    if (cutcount == 0)
        return false;

    cutq = LPcutList(head, cutcount);
    if (filter_primal)
        cutq.filter_primal(TG);

    return !cutq.empty();
}

}
}
