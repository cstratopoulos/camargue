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

Sep::SimpleDP::SimpleDP(CMR::Data::GraphGroup &graph_dat,
                        CMR::Data::KarpPartition &_kpart,
                        CMR::Data::BestGroup &best_dat,
                        CMR::Data::SupportGroup &supp_dat,
                        CMR::CutQueue<dominoparity> &_dp_q) try :
  candidates(graph_dat, best_dat, supp_dat), kpart(_kpart), dp_q(_dp_q)
  {} catch(const exception &e){
  cerr << e.what() << " constructing SimpleDP.\n";
  throw runtime_error("SimpleDP constructor failed.");
}

bool Sep::SimpleDP::find_cuts()
{
  runtime_error err("Problem in SimpleDP::find_cuts.");

  try {
    candidates.get_light_teeth();
    candidates.unmerged_weak_elim();
    candidates.merge_and_sort();
  } CMR_CATCH_PRINT_THROW("building and eliminating candidate teeth", err);

  for(int i = 0; i < kpart.num_parts(); ++i){
    try {
      CMR::DPwitness cutgraph(candidates, kpart[i]);
      cutgraph.simple_DP_sep(dp_q);
    } CMR_CATCH_PRINT_THROW("making a mini cutgraph sep call", err);
  }

  return(!dp_q.empty());
}

}
