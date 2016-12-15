#ifndef CMR_HYPERGRAPH_H
#define CMR_HYPERGRAPH_H

#include "cuts.hpp"
#include "cliq.hpp"

#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

class HyperGraph {
public:
    HyperGraph(CMR::Sep::CliqueBank &bank,
               const CCtsp_lpcut_in &cc_lpcut,
               const std::vector<int> &tour);
    
    ~HyperGraph();
    
private:
    char sense;
    double rhs;
    
    std::vector<CMR::Sep::Clique::Ptr> cliques;

    CMR::Sep::CliqueBank &source_bank;
};

class DominoCut {
public:
    DominoCut(CMR::Sep::CliqueBank &bank,
              CMR::dominoparity &dp_cut,
              const std::vector<int> &tour);

    ~DominoCut();

private:
    char sense;
    double rhs;

    CMR::Sep::Clique::Ptr handle;
    std::vector<std::pair<int, int>> nonneg_edges;
    std::vector<std::pair<int, CMR::Sep::Clique::Ptr>> teeth;

    CMR::Sep::CliqueBank &source_bank;
};

}
}

#endif
