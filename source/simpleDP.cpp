#include "simpleDP.hpp"
#include "util.hpp"
#include "err_util.hpp"
#include "config.hpp"

#include <iostream>
#include <memory>
#include <stdexcept>

using std::cout;
using std::cerr;

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
                        Sep::CutQueue<dominoparity> &_dp_q) try :
  candidates(active_tour, supp_dat), kpart(_kpart), dp_q(_dp_q)
{} catch (const exception &e) {
    cerr << e.what() << " constructing SimpleDP.\n";
    throw runtime_error("SimpleDP constructor failed.");
}

#ifndef CMR_USE_OMP
/////////////////////// SERIAL IMPLEMENTATION //////////////////////////////////
bool Sep::SimpleDP::find_cuts()
{

    runtime_error err("Problem in SimpleDP::find_cuts.");

    try {
        candidates.get_light_teeth();
        candidates.sort_by_root();
    } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);

    for (int i = 0; i < kpart.num_parts(); ++i) {
        CutQueue<dominoparity> mini_q;

        try {
            DPwitness cutgraph(candidates, kpart[i]);
            cutgraph.simple_DP_sep(mini_q);
        } CMR_CATCH_PRINT_THROW("making a mini cutgraph sep call", err);

        dp_q.splice(mini_q);
        // if(dp_q.size() >= 8)
        //     break;
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

    try {
        candidates.get_light_teeth();
        candidates.sort_by_root();
    } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);

    #pragma omp parallel for
    for (int i = 0; i < kpart.num_parts(); ++i) {
        if(at_capacity || caught_exception)
            continue;
        CutQueue<dominoparity> mini_q;

        try {
            DPwitness cutgraph(candidates, kpart[i]);

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
            dp_q.splice(mini_q);
            // if(dp_q.size() >= 25)
            //     at_capacity = true;
        }
    }

    if (caught_exception) {
        cerr << "OMP simple DP sep reported exception.\n";
        throw err;
    }

    return (!dp_q.empty());
}

#endif

}
