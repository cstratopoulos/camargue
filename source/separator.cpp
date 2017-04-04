#include "separator.hpp"
#include "cc_lpcuts.hpp"
#include "pool_sep.hpp"
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
                     Data::KarpPartition &kpart, int seed) try
    : core_edges(core_edges_), active_tour(active_tour_), supp_data(suppdata),
      karp_part(kpart),
      perm_elist(supp_data.support_elist), random_seed(seed)
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

    double gh2mt = 0.0;
    double f2mt = util::zeit();
    bool result = fast2m.find_cuts();
    f2mt = util::zeit() - f2mt;

    if (!result) {
        GHblossoms gh2m(perm_elist, supp_data.support_ecap, TG, fast2m_q);
        gh2mt = util::zeit();
        result = gh2m.find_cuts();
        gh2mt = util::zeit() - gh2mt;
    }

    if (verbose) {
        printf("\t%d fast blossoms, primal %d\t%.2fs\n",
               fast2m_q.size(), filter_primal, f2mt + gh2mt);
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
    if (supp_data.connected) {
        if (supp_data.in_subtour_poly()) {
            SimpleDP dominos(karp_part, active_tour, supp_data,
                             dp_q, random_seed);
            dominos.verbose = verbose;

            return dominos.find_cuts();
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

bool Separator::exsub_sep() try
{
    set_TG();
    ExactSub subtour(perm_elist, supp_data.support_ecap, TG, exsub_q);
    double exst = util::zeit();
    bool result = subtour.find_cuts();
    exst = util::zeit() - exst;

    if (verbose) {
        printf ("\t%d exact SECs\t%.2fs\n", exsub_q.size(), exst);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Separator::connect_sep failed.");
}

bool Separator::pool_sep(ExternalCuts &EC) try
{
    set_TG();

    PoolCuts pool_cuts(perm_elist, supp_data.support_ecap, TG, pool_q,
                       EC.cc_pool, random_seed);
    pool_cuts.filter_primal = filter_primal;

    double poolt = util::zeit();
    bool result = pool_cuts.find_cuts();
    poolt = util::zeit() - poolt;

    if (verbose) {
        printf("\t%d pool cuts\t%.2fs\n", pool_q.size(), poolt);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Separator::pool_sep failed");
}

bool Separator::tighten_pool(ExternalCuts &EC) try
{
    set_TG();

    PoolCuts pool_cuts(perm_elist, supp_data.support_ecap, TG, pool_q,
                       EC.cc_pool, random_seed);
    pool_cuts.filter_primal = filter_primal;

    double tightt = util::zeit();
    bool result = pool_cuts.tighten_pool();
    tightt = util::zeit() - tightt;

    if (result) {
        printf("\t%d tighten pool cuts\t%.2fs\n", pool_q.size(), tightt);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Separator::pool_sep failed");
}

bool Separator::consec1_sep(ExternalCuts &EC) try
{
    set_TG();

    PoolCuts pool_cuts(perm_elist, supp_data.support_ecap, TG, pool_q,
                       EC.cc_pool, random_seed);
    pool_cuts.filter_primal = filter_primal;

    double c1t = util::zeit();
    bool result = pool_cuts.find_consec1(EC.tightcuts);
    c1t = util::zeit() - c1t;

    if (result){
        printf("\t%d consec1 combs\t%.2fs\n", con1_q.size(), c1t);
        cout << flush;
    }

    return result;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Separator::consec1_sep failed");
}

bool Separator::local_sep(int chunk_sz, bool sphere) try
{
    set_TG();

    LocalCuts local_cuts(perm_elist, supp_data.support_ecap, TG, local_q,
                         random_seed);
    local_cuts.filter_primal = filter_primal;
    local_cuts.current_max = chunk_sz;
    local_cuts.spheres = sphere;

    double lct = util::zeit();
    bool result = local_cuts.find_cuts();
    lct = util::zeit() - lct;

    if (verbose) {
        if (!sphere) {
            printf("\t%d chunk %d local cuts\t%.2fs\n",
                   local_q.size(), chunk_sz, lct);
            cout << flush;
        } else {
            printf("\t%d chunk %d local cuts spheres\t%.2fs\n",
                   local_q.size(), chunk_sz, lct);
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
