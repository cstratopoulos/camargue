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
    HyperGraph(CliqueBank &bank,
               const CCtsp_lpcut_in &cc_lpcut,
               const std::vector<int> &tour);

    HyperGraph(CliqueBank &bank, ToothBank &tbank,
               const dominoparity &dp_cut, const double _rhs,
               const std::vector<int> &tour);

    /** Destruct a HyperGraph.
     * Goes through the source Cliquebank, decrementing reference count
     * on all Cliques used in the HyperGraph, purging Cliques from 
     * the CliqueBank if necessary. 
     */
    ~HyperGraph();

    enum Type : bool {
        Standard = true, Domino = false
    };

    Type cut_type() const { return static_cast<Type>(teeth.empty()); }
    
private:
    char sense;
    double rhs;
    
    std::vector<Clique::Ptr> cliques;
    std::vector<Tooth::Ptr> teeth;

    CliqueBank &source_bank;
    ToothBank *source_toothbank;
};


}
}

#endif
