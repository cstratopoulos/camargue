#include "pool_sep.hpp"
#include "err_util.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>

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
