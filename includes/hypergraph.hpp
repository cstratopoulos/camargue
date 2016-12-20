#ifndef CMR_HYPERGRAPH_H
#define CMR_HYPERGRAPH_H

#include "cut_structs.hpp"
#include "cliq.hpp"

#include <map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

/** Class for external representation of cuts added to the lp relaxation. */
class HyperGraph {
public:
    /** Construct a HyperGraph from a Concorde cut.
     * The cut corresponding to \p cc_lpcut will be represented using 
     * Clique pointers from \p bank, assuming that the cut was found with
     * \p tour as the resident best tour.
     */
    HyperGraph(CMR::Sep::CliqueBank &bank,
               const CCtsp_lpcut_in &cc_lpcut,
               const std::vector<int> &tour);

    /** Destruct a HyperGraph.
     * Goes through the source Cliquebank, decrementing reference count
     * on all Cliques used in the HyperGraph, purging Cliques from 
     * the CliqueBank if necessary. 
     */
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
              CMR::Sep::dominoparity &dp_cut, int rhs,
              const std::vector<int> &tour);

    ~DominoCut();

private:
    char sense;
    double rhs;

    CMR::Sep::Clique::Ptr handle;
    std::vector<CMR::Sep::Clique::Ptr> nonneg_edges;
    std::vector<std::pair<int, CMR::Sep::Clique::Ptr>> teeth;

    CMR::Sep::CliqueBank &source_bank;
};

/** Class for managing a list of cuts in an lp relaxation. */
class ExternalCuts {
    /** Construct ExternalCuts with a reference tour and perm. */
    ExternalCuts(const std::vector<int> &tour, const std::vector<int> &perm);

    /** Add a cut in Concorde format. */
    void add_cut(const CCtsp_lpcut_in &cc_lpcut,
                 const std::vector<int> &current_tour);

    /** Add a simple DP cut. */
    void add_cut(CMR::Sep::dominoparity &dp_cut,
                 const std::vector<int> &current_tour);

    /** Delete the cuts indicated by \p delset. 
     * If `delset[i] == -1` then `cuts[i]` will be deleted, else `cuts[i]` 
     * will become `cuts[delset[i]]`
     */
    void del_cuts(const std::vector<int> &delset,
                  int new_num_rows);

private:
    int next_row;
    
    CMR::Sep::CliqueBank clique_bank;
    
    std::map<int, CMR::Sep::HyperGraph> cuts;
    std::map<int, CMR::Sep::DominoCut> dp_cuts;
};

}
}

#endif
