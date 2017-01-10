#ifndef CMR_PRICE_UTIL_H
#define CMR_PRICE_UTIL_H

#include "util.hpp"

#include <array>

namespace CMR {

namespace Price {

constexpr int Nearest = 50; //<! Number of nearest edges to each node examined.

constexpr int AddBatch = 100; //<! Number added to CoreLP at a time.
constexpr int PoolSize = 1000; //<! Number of negative rc edges to keep in pool.
constexpr int EstBatch = 20000; //<! Max number of edges to estimate red cost.
constexpr int ScaleBatch = 3; //<! Scale factor for EstBatch.

constexpr double MaxPenalty = 0.10;

/// Return type for edge pricing routines.
enum class ScanStat {
    Partial, //!< Scanned some edges, found some with negative reduced cost.
    PartOpt, //!< Scanned some edges, found none with negative reduced cost.
    Full, //!< Scanned all edges, found some with negative reduced cost.
    FullOpt //!< Scanned all edges, none had negative reduced cost. 
};

struct PrEdge : EndPts {
    PrEdge() = default;
    
    PrEdge(int end0, int end1) : EndPts(end0, end1), redcost(1.0) {}

    PrEdge(int end0, int end1, double rc) : EndPts(end0, end1), redcost(rc) {}

    bool operator<(const PrEdge &rhs) const { return redcost < rhs.redcost; }
    
    double redcost;
};

}

}

#endif
