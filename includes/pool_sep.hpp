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

class PoolCuts : public CCsepBase {
public:
    PoolCuts(std::vector<int> &elist, std::vector<double> &ecap,
             TourGraph &TG, LPcutList &cutq, CCtsp_lpcuts *_pool, int seed)
        : CCsepBase(elist, ecap, TG, cutq), pool(_pool), random_seed(seed) {}

    bool find_cuts();

    bool tighten_pool();

    bool find_consec1(CCtsp_cuttree &tightcuts);

private:
    bool attempt_tighten();

    bool above_threshold(int num_paths);

    CCtsp_lpcuts *pool;
    int random_seed;
};

}
}

#endif
