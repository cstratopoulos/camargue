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

Separator::Separator(Data::GraphGroup &graphdata,
                     Data::BestGroup &bestdata,
                     Data::SupportGroup &suppdata,
                     Data::KarpPartition &kpart, Graph::TourGraph &_TG) try
    : max_total(8), running_total(0),
      graph_data(graphdata), best_data(bestdata), supp_data(suppdata),
      karp_part(kpart), TG(_TG), perm_elist(supp_data.support_elist)
{
    for (int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = best_data.perm[perm_elist[i]];
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator constructor failed.");
}

Separator::Separator(Data::GraphGroup &graphdata,
                     Data::BestGroup &bestdata,
                     Data::SupportGroup &suppdata,
                     Data::KarpPartition &kpart, Graph::TourGraph &_TG,
                     int round_limit) try
    : max_total(round_limit), running_total(0),
      graph_data(graphdata), best_data(bestdata), supp_data(suppdata),
      karp_part(kpart), TG(_TG), perm_elist(supp_data.support_elist)
{
    for (int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = best_data.perm[perm_elist[i]];
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator constructor failed.");
}

void ptr_reset(std::unique_ptr<Separator> &sep, Data::GraphGroup &g_dat,
               Data::BestGroup &b_dat, Data::SupportGroup &s_dat,
               Data::KarpPartition &kpart,
               Graph::TourGraph &TG)
{
    sep = util::make_unique<Separator>(g_dat, b_dat, s_dat, kpart, TG);
}

bool Separator::segment_sep() try
{
    SegmentCuts segments(perm_elist, supp_data.support_ecap, TG, seg_q);
    bool result = segments.find_cuts();

    running_total += seg_q.size();
    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::segment_sep failed.");
}

bool Separator::fast2m_sep() try
{
    FastBlossoms fast2m(perm_elist, supp_data.support_ecap, TG, fast2m_q);
    bool result = fast2m.find_cuts();

    running_total += fast2m_q.size();
    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::fast2m_sep failed.");
}

bool Separator::blkcomb_sep() try
{
    BlockCombs blkcomb(perm_elist, supp_data.support_ecap, TG, blkcomb_q);
    bool result = blkcomb.find_cuts();

    running_total += blkcomb_q.size();
    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::blkcomb_sep failed.");
}

bool Separator::exact2m_sep() try
{
    ExBlossoms ex2m(graph_data.core_graph.get_edges(),
                    best_data, supp_data,
                    ex2m_q);
    bool result = ex2m.find_cuts();

    running_total += ex2m_q.size();
    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::exact2m_sep failed.");
}

bool Separator::simpleDP_sep() try {
    if (supp_data.connected)
        if (supp_data.in_subtour_poly()) {
            SimpleDP dominos(graph_data, karp_part, best_data, supp_data,
                             dp_q);
            if (dominos.find_cuts()) {
                running_total += dp_q.size();
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

    running_total += connect_q.size();
    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::connect_sep failed.");
}


}
}
