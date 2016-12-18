/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief EXACT PRIMAL LIGHT SIMPLE DP SEPARATION
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "witness.hpp"
#include "tooth.hpp"
#include "datagroups.hpp"
#include "karp.hpp"
#include "process_cuts.hpp"

#include <vector>

namespace CMR {
namespace Sep {

/** Class for separation of simple domino parity inequalities. */
class SimpleDP {
public:
  SimpleDP(CMR::Data::GraphGroup &graph_dat,
           CMR::Data::KarpPartition &_kpart,
           CMR::Data::BestGroup &best_dat,
           CMR::Data::SupportGroup &supp_dat,
           CMR::Sep::CutQueue<dominoparity> &dp_q);

  bool find_cuts();

private:
  CMR::CandidateTeeth candidates;
  CMR::Data::KarpPartition &kpart;
  CMR::Sep::CutQueue<dominoparity> &dp_q;
};

}
}
