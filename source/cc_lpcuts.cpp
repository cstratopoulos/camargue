#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>
#include <utility>

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
    throw std::runtime_error("TourGraph constructor failed.");
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

///////////////////////////////////////////////////////////////////////////////
//                        CONCORDE SEPARATORS                                //
///////////////////////////////////////////////////////////////////////////////

bool SegmentCuts::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;

  if (CCtsp_segment_cuts(&head, &cutcount, TG.node_count(), ecap.size(),
                         &elist[0], &ecap[0]))
    throw std::runtime_error("CCtsp_segment_cuts failed.");

  cutq = LPcutList(head, cutcount);

  return(!cutq.empty());
}

bool ConnectCuts::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;

  if (CCtsp_connect_cuts(&head, &cutcount, TG.node_count(), ecap.size(),
			&elist[0], &ecap[0]))
    throw std::runtime_error("CCtsp_segment_cuts failed.");

  cutq = LPcutList(head, cutcount);

  return(!cutq.empty());
}

bool BlockCombs::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;

  if (CCtsp_block_combs(&head, &cutcount, TG.node_count(), ecap.size(),
		       &elist[0], &ecap[0], 1))
    throw std::runtime_error("CCtsp_block_combs failed.");

  if (cutcount == 0)
      return false;

  cutq = LPcutList(head, cutcount);

  if (filter_primal)
      cutq.filter_primal(TG);

  return(!cutq.empty());
}

bool FastBlossoms::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;

  if (CCtsp_fastblossom(&head, &cutcount, TG.node_count(), ecap.size(),
		       &elist[0], &ecap[0]))
    throw std::runtime_error("CCtsp_fastblossom failed.");

  if (cutcount == 0) {
      if (CCtsp_ghfastblossom(&head, &cutcount, TG.node_count(), ecap.size(),
                              &elist[0], &ecap[0]))
          throw std::runtime_error("CCtsp_ghfastblossom failed.");
  }

  cutq = LPcutList(head, cutcount);
  if (filter_primal)
      cutq.filter_primal(TG);

  if (!cutq.empty())
      return true;

  return(!cutq.empty());
}

}
}
