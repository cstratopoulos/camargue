/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief EXACT PRIMAL LIGHT SIMPLE DP SEPARATION
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "DPgraph.hpp"
#include "tooth.hpp"
#include "datagroups.hpp"
#include "cuts.hpp"

namespace CMR {

template<> class Cut<dominoparity> {
public:
  Cut<dominoparity>(CMR::Data::GraphGroup &graph_dat,
		    CMR::Data::BestGroup &best_dat,
		    CMR::Data::SupportGroup &supp_dat,
		    CMR::CutQueue<dominoparity> &_dp_q) :
    candidates(graph_dat, best_dat, supp_dat),
    dp_q(_dp_q) {}

  int cutcall();

protected:
  int separate();
  int add_cuts();

private:
  CMR::CandidateTeeth candidates;
  CMR::CutQueue<CMR::dominoparity> &dp_q;
  //CMR::CutQueue<CMR::externalDP> &external_q;
};

}
