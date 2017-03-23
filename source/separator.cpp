#include "separator.hpp"
#include "cc_lpcuts.hpp"
#include "simpleDP.hpp"
#include "blossoms.hpp"
#include "err_util.hpp"

#include <vector>
#include <stdexcept>
#include <iostream>

#include <cstdio>

using std::vector;

using std::exception;
using std::runtime_error;
using std::logic_error;

using std::cout;
using std::cerr;
using std::endl;
using std::flush;


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

    double st = util::zeit();
    bool result = segments.find_cuts();
    st = util::zeit() - st;

    if (verbose) {
        printf("\t%d segment cuts\t%.2fs\n", seg_q.size(), st);
        cout << flush;
    }

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

    double f2mt = util::zeit();
    bool result = fast2m.find_cuts();
    f2mt = util::zeit() - f2mt;

    if (verbose) {
        printf("\t%d fast blossoms, primal %d\t%.2fs\n",
               fast2m_q.size(), filter_primal, f2mt);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::fast2m_sep failed.");
}

bool Separator::blkcomb_sep() try
{
    set_TG();
    BlockCombs blkcomb(perm_elist, supp_data.support_ecap, TG, blkcomb_q);

    double blkt = util::zeit();
    bool result = blkcomb.find_cuts();
    blkt = util::zeit() - blkt;

    if (verbose) {
        printf("\t%d block combs, primal %d\t%.2fs\n",
               blkcomb_q.size(), filter_primal, blkt);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Separator::blkcomb_sep failed.");
}

bool Separator::exact2m_sep() try
{
    set_TG();
    ExBlossoms ex2m(core_edges, active_tour, supp_data, ex2m_q);

    Timer e2mt("Primal blossoms");
    e2mt.start();
    bool result = ex2m.find_cuts();
    e2mt.stop();

    if (verbose) {
        cout << "\t" << ex2m_q.size() << " primal blossoms" << endl;
#ifdef CMR_USE_OMP
        bool cpu = true;
#else
        bool cpu = false;
#endif
        e2mt.report(cpu);
    }

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
            dominos.verbose = verbose;

            return dominos.find_cuts();
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
    double cont = util::zeit();
    bool result = subtour.find_cuts();
    cont = util::zeit() - cont;

    if (verbose) {
        printf ("\t%d connect cuts\t%.2fs\n", connect_q.size(), cont);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Separator::connect_sep failed.");
}

bool Separator::local_sep() try
{
    set_TG();

    LocalCuts local_cuts(perm_elist, supp_data.support_ecap, TG, local_q);
    local_cuts.current_max = lc_chunk;
    local_cuts.spheres = lc_sphere;

    double lct = util::zeit();
    bool result = local_cuts.find_cuts();
    lct = util::zeit() - lct;

    if (verbose) {
        if (!lc_sphere) {
            printf("\t%d chunk %d local cuts\t%.2fs\n",
                   local_q.size(), lc_chunk, lct);
            cout << flush;
        } else {
            printf("\t%d chunk %d local cuts spheres\t%.2fs\n",
                   local_q.size(), lc_chunk, lct);
            cout << flush;
        }
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Separator::local_sep failed.");
}


}
}
