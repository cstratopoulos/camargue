#include "separator.hpp"
#include "cc_lpcuts.hpp"
#include "simpleDP.hpp"
#include "blossoms.hpp"
#include "err_util.hpp"

#include <vector>
#include <stdexcept>
#include <iostream>

using std::vector;

using std::exception;
using std::runtime_error;
using std::logic_error;

using std::cout;
using std::cerr;


namespace CMR {
namespace Sep {

Separator::Separator(const vector<Graph::Edge> &core_edges_,
                     const LP::ActiveTour &active_tour_,
                     Data::SupportGroup &suppdata,
                     Data::KarpPartition &kpart) try
    : core_edges(core_edges_), active_tour(active_tour_), supp_data(suppdata),
      karp_part(kpart),
      perm_elist(supp_data.support_elist)
{
    for (int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = active_tour.tour_perm()[perm_elist[i]];
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator constructor failed.");
}

/**
 * @warning This method must be called before constructing any any separator
 * that derives from ConcordeSeparator, i.e., SegmentCuts, FastBlossoms,
 * BlockCombs, and ConnectCuts.
 */
void Separator::set_TG()
{
    TG = TourGraph(active_tour.edges(), core_edges, active_tour.tour_perm());
}

bool Separator::segment_sep() try
{
    set_TG();
    SegmentCuts segments(perm_elist, supp_data.support_ecap, TG, seg_q);
    bool result = segments.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::segment_sep failed.");
}

bool Separator::fast2m_sep() try
{
    set_TG();
    FastBlossoms fast2m(perm_elist, supp_data.support_ecap, TG, fast2m_q);

    fast2m.filter_primal = filter_primal;

    bool result = fast2m.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::fast2m_sep failed.");
}

bool Separator::blkcomb_sep() try
{
    set_TG();
    BlockCombs blkcomb(perm_elist, supp_data.support_ecap, TG, blkcomb_q);
    bool result = blkcomb.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::blkcomb_sep failed.");
}

bool Separator::exact2m_sep() try
{
    set_TG();
    ExBlossoms ex2m(core_edges, active_tour, supp_data, ex2m_q);

    bool result = ex2m.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::exact2m_sep failed.");
}

bool Separator::simpleDP_sep() try {
    if (supp_data.connected)
        if (supp_data.in_subtour_poly()) {
            SimpleDP dominos(karp_part, active_tour, supp_data,
                             dp_q);
            if (dominos.find_cuts()) {
                return true;
            }
        }
    return false;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::simpleDP_sep failed.");
}

bool Separator::connect_sep() try
{
    set_TG();
    ConnectCuts subtour(perm_elist, supp_data.support_ecap, TG, connect_q);
    bool result = subtour.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::connect_sep failed.");
}


}
}
