#include "pool_sep.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cstdio>

extern "C" {
#include <concorde/INCLUDE/consec1.h>
#include <concorde/INCLUDE/cuttree.h>
}

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::exception;

using std::vector;

namespace CMR {

namespace Eps = Epsilon;

namespace Sep {

bool PoolCuts::find_cuts()
{
    int cutcount = 0;
    double max_viol = 0.0;
    CCtsp_lpcut_in *head = NULL;
    CCrandstate rstate;

    CCutil_sprand(random_seed, &rstate);

    if (CCtsp_search_cutpool(pool, &head, &cutcount, &max_viol,
                             TG.node_count(), ecap.size(), &elist[0],
                             &ecap[0], 0, &rstate))
        throw runtime_error("CCtsp_search_cutpool failed");

    if (cutcount == 0)
        return false;

    cutq = LPcutList(head, cutcount);
    if (filter_primal)
        cutq.filter_primal(TG);

    return !cutq.empty();
}

bool PoolCuts::tighten_pool()
{
    runtime_error err("Problem in PoolCuts::tighten_pool");

    try {
        if (!attempt_tighten())
            return false;
    } CMR_CATCH_PRINT_THROW("calling attempt_sep", err);

    int cutcount = 0;
    double max_viol = 0.0;
    CCtsp_lpcut_in *head = NULL;
    CCrandstate rstate;
    CCtsp_tighten_info _ti_unused;

    CCtsp_init_tighten_info(&_ti_unused);
    CCutil_sprand(random_seed, &rstate);

    if (CCtsp_tighten_lp(pool, &_ti_unused, &head, &cutcount,
                         TG.node_count(), ecap.size(),
                         &elist[0], &ecap[0], 0.1, 250,
                         &max_viol, &rstate)) {
        cerr << "CCtsp_tighten_lp failed" << endl;
        throw err;
    }

    if (cutcount == 0)
        return false;

    cutq = LPcutList(head, cutcount);
    if (filter_primal)
        cutq.filter_primal(TG);

    return !cutq.empty();
}

bool PoolCuts::find_consec1(CCtsp_cuttree &tight_cuts)
{
    runtime_error err("Problem in PoolCuts::find_consec1");

    int ecount = ecap.size();

    if (CCpq_cuttree_improve_quick(&tight_cuts, pool, ecount, &elist[0],
                                   &ecap[0])) {
        cerr << "CCpq_cuttree_improve_quick failed" << endl;
        throw err;
    }

    int cutcount = 0;
    CCtsp_lpcut_in *head = NULL;

    if (CCpq_consecutiveones(&head, &cutcount, &tight_cuts, pool, ecount,
                             &elist[0], &ecap[0])) {
        cerr << "CCpq_consecutiveones failed" << endl;
        throw err;
    }

    if (cutcount == 0)
        return false;

    cutq = LPcutList(head, cutcount);
    if (filter_primal)
        cutq.filter_primal(TG);

    return !cutq.empty();
}

bool PoolCuts::find_tour_tight()
{
    runtime_error err("Problem in PoolCuts::find_tour_tight");

    int ncount = TG.node_count();
    CCtsp_lpgraph *tour_lpg = TG.pass_ptr();
    int ecount = tour_lpg->ecount;

    vector<int> tour_elist;

    try { tour_elist.reserve(ecount); }
    CMR_CATCH_PRINT_THROW("reserving tour elist", err);

    for (int i = 0; i < ecount; ++i) {
        tour_elist.push_back(tour_lpg->edges[i].ends[0]);
        tour_elist.push_back(tour_lpg->edges[i].ends[1]);
    }

    vector<double> cutvals;

    try { cutvals.resize(pool->cutcount); }
    CMR_CATCH_PRINT_THROW("allocating cutvals", err);

    if (CCtsp_price_cuts(pool, ncount, ecount, &tour_elist[0],
                         TG.tour_array(), &cutvals[0]))
        throw runtime_error("CCtsp_search_cutpool failed");

    int num_tight = std::count(cutvals.begin(), cutvals.end(), 0.0);
    if (num_tight == 0) {
        cout << "No tight cuts found : (" << endl;
        return false;
    } else {
        cout << num_tight << " cuts are tight at tour" << endl;
    }

    int num_added = 0;

    for (int i = 0; i < pool->cutcount && num_added <= 500; ++i)
        if (cutvals[i] == 0.0) {
            CCtsp_lpcut_in *newc = CC_SAFE_MALLOC(1, CCtsp_lpcut_in);
            if (newc == NULL) {
                cerr << "Out of memory for new cut" << endl;
                throw err;
            }

            if (CCtsp_lpcut_to_lpcut_in(pool, &pool->cuts[i], newc)) {
                cerr << "CCtsp_lpcut_to_lpcut_in failed" << endl;
                CC_FREE(newc, CCtsp_lpcut_in);
                throw err;
            }

            cutq.push_front(newc);
            ++num_added;
        }

    return true;
}

bool PoolCuts::above_threshold(int num_paths)
{
    return num_paths > 0.2 * TG.node_count();
}

bool PoolCuts::attempt_tighten()
{
    CC_SRKgraph G;
    auto cleanup = util::make_guard([&G]{ CCcut_SRK_free_graph(&G); });

    int num_paths = 0;
    int ncount = TG.node_count();

    CCcut_SRK_init_graph(&G);
    if (CCcut_SRK_buildgraph(&G, ncount, ecap.size(),
                             &elist[0], &ecap[0]))
        throw runtime_error("CCcut_SRK_buildgraph failed");

    CCcut_SRK_increment_marker(&G);

    if (CCcut_SRK_defluff(&G))
        throw runtime_error("CCcut_SRK_defluff failed");

    CCcut_SRK_identify_paths_to_edges(&G, &num_paths, 0);

    return above_threshold(num_paths);
}

}
}
