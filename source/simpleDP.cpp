#include "simpleDP.hpp"
#include "util.hpp"
#include "err_util.hpp"
#include "timer.hpp"
#include "config.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;

using std::unique_ptr;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

/**
 * Constructs a separator from a core LP graph and LP solution.
 * @param[in] graph_dat the info of the core edge set/graph.
 * @param[in] _kpart the KarpPartition for separating over partitioned
 * DPwitness graphs.
 * @param[in] supp_dat the SupportGroup for the current LP solution.
 * @param[in] _dp_q the CutQueue where all found cuts will be stored.
 */
Sep::SimpleDP::SimpleDP(Data::KarpPartition &_kpart,
                        const LP::ActiveTour &active_tour,
                        Data::SupportGroup &supp_dat,
                        Sep::CutQueue<dominoparity> &_dp_q,
                        int seed) try :
    candidates(active_tour, supp_dat), kpart(_kpart), dp_q(_dp_q),
    random_seed(seed)
    {} catch (const exception &e) {
    cerr << e.what() << " constructing SimpleDP.\n";
    throw runtime_error("SimpleDP constructor failed.");
}

#ifndef CMR_USE_OMP
/////////////////////// SERIAL IMPLEMENTATION //////////////////////////////////
bool Sep::SimpleDP::find_cuts()
{

    runtime_error err("Problem in SimpleDP::find_cuts.");

    Timer find_total("Serial simple DP sep");
    Timer find_cands("Finding cand teeth", &find_total);

    find_total.start();
    find_cands.start();

    if (verbose)
        cout << "\tSerial simple DP sep....\n";

    try {
        candidates.get_light_teeth();
        candidates.sort_by_root();
        if (verbose) {
            int numfound = 0;
            for (const auto &v : candidates.light_teeth)
                numfound += v.size();
            cout << "\t\tGot " << numfound
                 << " light teeth" << endl;
        }

    } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);

    find_cands.stop();

    Timer search_wit("Make/search witness", &find_total);
    search_wit.start();
    for (int i = 0; i < kpart.num_parts(); ++i) {
        CutQueue<dominoparity> mini_q;


        try {
            DPwitness cutgraph(candidates, kpart[i], random_seed);
            cutgraph.simple_DP_sep(mini_q);
        } CMR_CATCH_PRINT_THROW("making a mini cutgraph sep call", err);

        if (verbose)
            cout << "\t\t" << mini_q.size() << " cuts from partition "
                 << i << "\n";
        dp_q.splice(mini_q);
        if(dp_q.size() >= 250) {
            cout << "\t\tBreaking early.\n";
            break;
        }
    }

    search_wit.stop();
    find_total.stop();

    if (verbose) {
        cout << "\t" << dp_q.size() << " total simple DP cuts enqueued."
             << endl;
        if (verbose > 1) {
            candidates.profile();
            find_cands.report(false);
            search_wit.report(false);
            find_total.report(false);
        }
    }

    return(!dp_q.empty());
}

#else
/////////////////////// OMP PARALLEL IMPLEMENTATION ////////////////////////////
bool Sep::SimpleDP::find_cuts()
{
    runtime_error err("Problem in SimpleDP::find_cuts.");

    bool at_capacity = false;
    bool caught_exception = false;

    Timer find_total("Parallel simple DP sep");
    Timer find_cands("finding candid teeth", &find_total);

    find_total.start();
    find_cands.start();
    try {
        candidates.get_light_teeth();
        candidates.sort_by_root();
        if (verbose) {
            int numfound = 0;
            for (const auto &v : candidates.light_teeth)
                numfound += v.size();
            cout << "\t\tGot " << numfound
                 << " light teeth" << endl;
        }
    } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);
    find_cands.stop();

    if (verbose)
        cout << "Parallel search over witness graphs" << endl;

    Timer search_wit("make/search witness", &find_total);
    search_wit.start();
    #pragma omp parallel for
    for (int i = 0; i < kpart.num_parts(); ++i) {
        if(at_capacity || caught_exception)
            continue;
        CutQueue<dominoparity> mini_q;

        try {
            DPwitness cutgraph(candidates, kpart[i], random_seed);

            cutgraph.simple_DP_sep(mini_q);
        } catch (const exception &e) {
            #pragma omp critical
            {
                cerr << "Caught " << e.what() << " in witness subproblem.\n";
                caught_exception = true;
            }
        }

        if(caught_exception)
            continue;

        #pragma omp critical
        {
        if (verbose)
            cout << "\t" << mini_q.size() << " cuts from partition "
                 << i << "\n";

            dp_q.splice(mini_q);

            if(dp_q.size() >= 250 && !at_capacity) {
                if (verbose)
                    cout << "DP q has size " << dp_q.size() << ", "
                         << "terminating on part number "
                         << i << endl;
                at_capacity = true;
            }
        }
    }

    search_wit.stop();
    find_total.stop();

    if (caught_exception) {
        cerr << "OMP simple DP sep reported exception.\n";
        throw err;
    }

    if (verbose) {
        cout << dp_q.size() << " total simple DP cuts enqueued." << endl;
        if (verbose > 1) {
            candidates.profile();
            find_cands.report(true);
            search_wit.report(true);
            find_total.report(true);
        }
    }

    return (!dp_q.empty());
}

#endif

}
