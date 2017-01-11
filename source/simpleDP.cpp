#include "simpleDP.hpp"
#include "util.hpp"
#include "err_util.hpp"
#include "tests.hpp"

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

Sep::SimpleDP::SimpleDP(Data::GraphGroup &graph_dat,
                        Data::KarpPartition &_kpart,
                        Data::BestGroup &best_dat,
                        Data::SupportGroup &supp_dat,
                        Sep::CutQueue<dominoparity> &_dp_q) try :
  candidates(graph_dat, best_dat, supp_dat), kpart(_kpart), dp_q(_dp_q)
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
    candidates.unmerged_weak_elim();
    candidates.merge_and_sort();
  } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);

  for (int i = 0; i < kpart.num_parts(); ++i) {
      CutQueue<dominoparity> mini_q(25);
      
      try {
          DPwitness cutgraph(candidates, kpart[i]);      
          cutgraph.simple_DP_sep(mini_q);
      } CMR_CATCH_PRINT_THROW("making a mini cutgraph sep call", err);

      dp_q.splice(mini_q);
      if(dp_q.size() >= 8)
          break;
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
        candidates.unmerged_weak_elim();
        candidates.merge_and_sort();
    } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);

    #pragma omp parallel for
    for (int i = 0; i < kpart.num_parts(); ++i) {
        if(at_capacity || caught_exception)
            continue;
        CutQueue<dominoparity> mini_q(25);
        
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
            if(dp_q.size() >= 25)
                at_capacity = true;
        }
    }

    return (!dp_q.empty());
}

#endif

}
