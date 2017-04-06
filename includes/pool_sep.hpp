/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Separating over a pool of cuts.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_POOL_SEP_H
#define CMR_POOL_SEP_H

#include "cc_lpcuts.hpp"
#include "hypergraph.hpp"

#include <unordered_map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

/// Cut pool separation.
class PoolCuts : public CCsepBase {
public:
    PoolCuts(std::vector<int> &elist, std::vector<double> &ecap,
             TourGraph &TG, LPcutList &cutq, CCtsp_lpcuts *_pool, int seed)
        : CCsepBase(elist, ecap, TG, cutq), pool(_pool), random_seed(seed) {}

    /// Search the pool for violated cuts.
    bool find_cuts();

    /// Try to obtain violated cuts by tightening cuts in the pool.
    bool tighten_pool();

    /// Find combs from consecutive ones.
    bool find_consec1(CCtsp_cuttree &tightcuts);

    /// Search the pool just for cuts that are tight at the current tour.
    bool find_tour_tight();

private:
    /// Based on the ecap vector, should we attempt tightening.
    bool attempt_tighten();

    /// Threshold used by PoolCuts#attempt_tighten.
    bool above_threshold(int num_paths);

    CCtsp_lpcuts *pool;
    int random_seed;
};

}
}

#endif
