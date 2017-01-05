#ifndef CMR_PRICE_UTIL_H
#define CMR_PRICE_UTIL_H

namespace CMR {

namespace Price {

constexpr int BatchSize = 100; //<! Number of edges added to CoreLP at a time.
constexpr int QueueMax = 20000; //<! Max number of edges in queue for addition.

/// Return type for edge pricing routines.
enum class ScanStat {
    Partial, //!< Scanned some edges, found some with negative reduced cost.
    PartOpt, //!< Scanned some edges, found none with negative reduced cost.
    Full, //!< Scanned all edges, found some with negative reduced cost.
    FullOpt //!< Scanned all edges, none had negative reduced cost. 
};

struct edge {
    edge() = default;
    
    edge(int end0, int end1) :
        end{end0 < end1 ? end0 : end1, end1 > end0 ? end1 : end0},
        redcost(1.0) {}
    
    edge(int end0, int end1, double rc) :
        end{end0 < end1 ? end0 : end1, end1 > end0 ? end1 : end0},
        redcost(rc) {}

    bool operator<(const edge &rhs) const { return redcost > rhs.redcost; }
    
    std::array<int, 2> end;
    double redcost;
};

}

}

#endif
