#ifndef CMR_HYPERGRAPH_H
#define CMR_HYPERGRAPH_H

#include "cut_structs.hpp"
#include "cliq.hpp"

#include <map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

/** Class for external representation of a cut added to the lp relaxation. */
class HyperGraph {
public:
    HyperGraph() = default;
    
    /** Construct a HyperGraph from a Concorde cut.
     * The cut corresponding to \p cc_lpcut will be represented using 
     * Clique pointers from \p bank, assuming that the cut was found with
     * \p tour as the resident best tour.
     */
    HyperGraph(CliqueBank &bank,
               const CCtsp_lpcut_in &cc_lpcut,
               const std::vector<int> &tour);

    /** Construct a HyperGraph from a simple tooth inequality. */
    HyperGraph(CliqueBank &bank, ToothBank &tbank,
               const dominoparity &dp_cut, const double _rhs,
               const std::vector<int> &tour);

    ~HyperGraph(); /**< Destruct and decrement/delete Clique/Tooth refs. */

    enum Type : bool {
        Standard = true, Domino = false
    };

    Type cut_type() const { return static_cast<Type>(teeth.empty()); }

    /** Get the coefficient of an edge specified by endpoints. */
    double get_coeff(int end0, int end1) const;

    friend class ExternalCuts;
    
private:
    char sense;
    double rhs;
    
    std::vector<Clique::Ptr> cliques;
    std::vector<Tooth::Ptr> teeth;

    CliqueBank *source_bank;
    ToothBank *source_toothbank;
};

/** The external storage of a collection of HyperGraph cuts in a Relaxation. */
class ExternalCuts {
public:
    /** Construct ExternalCuts with reference tour and perm. */
    ExternalCuts(const std::vector<int> &tour, const std::vector<int> &perm);

    /** Add a Concorde cut. */
    void add_cut(const CCtsp_lpcut_in &cc_lpcut,
                 const std::vector<int> &current_tour);

    /** Add a simple DP cut. */
    void add_cut(const dominoparity &dp_cut, const double rhs,
                 const std::vector<int> &current_tour);

    /** Delete the cuts specified by \p delset.
     * Erase each entry of `cuts[i]` for which `delset[i + node_count] == -1`.
     */
    void del_cuts(const std::vector<int> &delset);

    /** Return a cut corresponding to a row number index from the lp. */
    const HyperGraph &get_cut(int lp_rownum) const {
        return cuts[lp_rownum - node_count];
    }

private:
    const int node_count;

    CliqueBank clique_bank;
    ToothBank tooth_bank;

    std::vector<HyperGraph> cuts;
};

}
}

#endif
