#ifndef CMR_PRICE_UTIL_H
#define CMR_PRICE_UTIL_H

namespace CMR {

namespace Price {

constexpr int Nearest = 50; //<! Number of nearest edges to each node examined.

constexpr int AddBatch = 100; //<! Number added to CoreLP at a time.
constexpr int PoolSize = 1000; //<! Number of negative rc edges to keep in pool.
constexpr int EstBatch = 20000; //<! Max number of edges to estimate red cost.
constexpr int ScaleBatch = 3; //<! Scale factor for EstBatch. 

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
        redcost(1.0)
  {
    end[0] = end0 < end1 ? end0 : end1;
    end[1] = end1 > end0 ? end1 : end0;
  }
    
    edge(int end0, int end1, double rc) :
        redcost(rc)
  {
    end[0] = end0 < end1 ? end0 : end1;
    end[1] = end1 > end0 ? end1 : end0;
  }

    bool operator<(const edge &rhs) const { return redcost > rhs.redcost; }
    
    std::array<int, 2> end;
    double redcost;
};

}

}

#endif
