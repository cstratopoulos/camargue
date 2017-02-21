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
                     Data::BestGroup &bestdata,
                     Data::SupportGroup &suppdata,
                     Data::KarpPartition &kpart) try
    : core_edges(core_edges_), best_data(bestdata), supp_data(suppdata),
      karp_part(kpart),
      TG(bestdata.best_tour_edges, core_edges_, bestdata.perm),
      perm_elist(supp_data.support_elist)
{
    for (int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = best_data.perm[perm_elist[i]];
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator constructor failed.");
}

bool Separator::segment_sep() try
{
    SegmentCuts segments(perm_elist, supp_data.support_ecap, TG, seg_q);
    bool result = segments.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::segment_sep failed.");
}

bool Separator::fast2m_sep() try
{
    FastBlossoms fast2m(perm_elist, supp_data.support_ecap, TG, fast2m_q);
    bool result = fast2m.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::fast2m_sep failed.");
}

bool Separator::blkcomb_sep() try
{
    BlockCombs blkcomb(perm_elist, supp_data.support_ecap, TG, blkcomb_q);
    bool result = blkcomb.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::blkcomb_sep failed.");
}

bool Separator::exact2m_sep() try
{
    ExBlossoms ex2m(core_edges, best_data, supp_data, ex2m_q);
    bool result = ex2m.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::exact2m_sep failed.");
}

bool Separator::simpleDP_sep() try {
    if (supp_data.connected)
        if (supp_data.in_subtour_poly()) {
            SimpleDP dominos(karp_part, best_data, supp_data,
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

bool Separator::connect_sep() try {
    ConnectCuts subtour(perm_elist, supp_data.support_ecap, TG, connect_q);
    bool result = subtour.find_cuts();

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::connect_sep failed.");
}


}
}
