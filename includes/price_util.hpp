#ifndef CMR_PRICE_UTIL_H
#define CMR_PRICE_UTIL_H

#include "util.hpp"

#include <array>
#include <iostream>
#include <string>

namespace CMR {

namespace Price {

constexpr int Nearest = 50; //<! Number of nearest edges to each node examined.

constexpr int AddBatch = 100; //<! Number added to CoreLP at a time.
constexpr int PoolSize = 1000; //<! Number of negative rc edges to keep in pool.
constexpr int EstBatch = 20000; //<! Max number of edges to estimate red cost.
constexpr int ScaleBatch = 3; //<! Scale factor for EstBatch.

constexpr int f64Batch = 5000000; //!< Max number to generate during exact lb.

constexpr double MaxPenalty = 0.10;

/// Return type for edge pricing routines.
enum class ScanStat {
    Partial, //!< Scanned some edges, found some with negative reduced cost.
    PartOpt, //!< Scanned some edges, found none with negative reduced cost.
    Full, //!< Scanned all edges, found some with negative reduced cost.
    FullOpt //!< Scanned all edges, none had negative reduced cost. 
};


inline std::ostream &operator<<(std::ostream &os, ScanStat stat)
{
    std::string out;
    switch (stat) {
    case ScanStat::Partial:
        out = "Found edges in Partial";
        break;
    case ScanStat::PartOpt:
        out = "Optimal for Partial";
        break;
    case ScanStat::Full:
        out = "Found edges in Full";
        break;
    case ScanStat::FullOpt:
        out = "Optimal for Full";
        break;
    }

    os << out;
    return os;
}


template <typename numtype>
struct PrEdge : EndPts {
    PrEdge() = default;
    
    PrEdge(int end0, int end1) : EndPts(end0, end1), redcost(1.0) {}

    PrEdge(int end0, int end1, numtype rc) : EndPts(end0, end1), redcost(rc) {}

    bool operator<(const PrEdge &rhs) const { return redcost < rhs.redcost; }
    
    numtype redcost;
};

}

}

#endif
