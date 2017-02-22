/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Primal light simple DP separation.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "active_tour.hpp"
#include "witness.hpp"
#include "tooth.hpp"
#include "datagroups.hpp"
#include "karp.hpp"
#include "process_cuts.hpp"

#include <vector>

namespace CMR {
namespace Sep {

/// Separating primal simple domino parity inequalities.
class SimpleDP {
public:
    /// Construct a separator to separate over partitioned DPwitness graphs.
    SimpleDP(Data::KarpPartition &_kpart,
             const LP::ActiveTour &active_tour,
             Data::SupportGroup &supp_dat,
             Sep::CutQueue<dominoparity> &_dp_q);

    bool find_cuts(); //!< Separator invocation, returns true iff cuts found.

private:
    CandidateTeeth candidates;
    Data::KarpPartition &kpart;
    CutQueue<dominoparity> &dp_q;
};

}
}
