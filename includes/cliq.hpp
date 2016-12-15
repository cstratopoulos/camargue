#ifndef CMR_CLIQ_H
#define CMR_CLIQ_H

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace CMR {
namespace Sep {

struct segment {
    segment() = default;
    segment(int lo, int hi) : start(lo), end(hi) {}

    int size() const { return end - start + 1; }

    bool operator>(const segment &rhs) const
    {
        return std::make_tuple(size(), start, end) >
            std::make_tuple(rhs.size(), rhs.start, rhs.end);
    }

    bool operator==(const segment &rhs) const
    { return (start == rhs.start) && (end == rhs.end); }
    
    int start;
    int end;
};

/** Class for storing segment lists representing edges of a hypergraph. */
class Clique {
public:
    Clique() = default;
    
    Clique(const CCtsp_lpclique &cc_cliq,
           const std::vector<int> &saved_tour,
           const std::vector<int> &saved_perm,
           const std::vector<int> & current_tour);

    Clique(std::vector<int> &nodes,
           const std::vector<int> &perm);

    using Ptr = std::shared_ptr<Clique>;

    int seg_count() const { return seglist.size(); }

    const std::vector<CMR::Sep::segment> &seg_list() const { return seglist; }
    std::vector<int> node_list(const std::vector<int> &saved_tour) const;
    

    bool operator==(const Clique &rhs) const { return seglist == rhs.seglist; }
    
private:
    std::vector<CMR::Sep::segment> seglist;
};

}
}

namespace std {

/** Partial specialization of std::hash taken from CCtsp_hashclique. */
template<>
struct hash<CMR::Sep::Clique> {
    size_t operator()(const CMR::Sep::Clique &clq) const
    {
        size_t val = 0;
        
        for (const CMR::Sep::segment &seg : clq.seg_list())
            val = (val * 65537) + (seg.start * 4099) + seg.end;

        return val;
    }
};
}

namespace CMR {
namespace Sep {

/** Storage of a repository of Cliques, for use in building a HyperGraph. */
class CliqueBank {
public:
    CliqueBank(const std::vector<int> &tour, const std::vector<int> &perm);

    CMR::Sep::Clique::Ptr add_clique(const CMR::Sep::Clique &clq);
    
    CMR::Sep::Clique::Ptr add_clique(const CCtsp_lpclique &cc_cliq,
                                      const std::vector<int> &tour);
    
    CMR::Sep::Clique::Ptr add_clique(std::vector<int> &nodes);
    
    void del_clique(CMR::Sep::Clique::Ptr &clq_ptr);

    using Itr = std::unordered_map<CMR::Sep::Clique,
                                    CMR::Sep::Clique::Ptr>::iterator;

    int size() const { return bank.size(); }
    
private:
    const std::vector<int> saved_tour;
    const std::vector<int> saved_perm;
    std::unordered_map<CMR::Sep::Clique, CMR::Sep::Clique::Ptr> bank;
};


}
}

#endif
