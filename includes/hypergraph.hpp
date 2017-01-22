#ifndef CMR_HYPERGRAPH_H
#define CMR_HYPERGRAPH_H

#include "cut_structs.hpp"
#include "lp_interface.hpp"
#include "cliq.hpp"
#include "price_util.hpp"

#include <map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

/// External representation of a cut added to the lp relaxation.
class HyperGraph {
public:
    HyperGraph() = default;
    

    /// Construct a HyperGraph from a Concorde cut.
    HyperGraph(CliqueBank &bank,
               const CCtsp_lpcut_in &cc_lpcut,
               const std::vector<int> &tour);

    /// Construct a HyperGraph from a simple tooth inequality.
    HyperGraph(CliqueBank &bank, ToothBank &tbank,
               const dominoparity &dp_cut, const double _rhs,
               const std::vector<int> &tour);

    ~HyperGraph(); //!< Destruct and decrement/delete Clique/Tooth refs.

    /// Enumeration for the types of HyperGraph inequalities.
    enum Type {
        Domino = 0, //<! A domino parity inequality.
        Subtour = 1, //<! An SEC.
        Comb = 2 //<! A comb-like constraint.
    };

    Type cut_type() const; //!< Find the Type of this cut.

    /// Get the coefficient of an edge specified by endpoints.
    double get_coeff(int end0, int end1) const;

    /// Get sparse coefficient row for a list of edges.
    void get_coeffs(const std::vector<Price::PrEdge> &edges,
                    std::vector<int> &rmatind,
                    std::vector<double> &rmatval) const;

    char get_sense() const { return sense; } //!< Get the sense of the cut.    
    double get_rhs() const { return rhs; } //!< Get the rhs of the cut.

    
    const std::vector<Clique::Ptr> &get_cliques() const
        { return cliques; } //!< Constant ref to the vector of Clique refs.

    const std::vector<Tooth::Ptr> &get_teeth() const
        { return teeth; } //!< Constant ref to the vector of Tooth refs.

    friend class ExternalCuts;
    
private:
    char sense; //<! The inequality sense of the cut.
    double rhs; //<! The righthand-side of the cut.
    
    std::vector<Clique::Ptr> cliques; //<! The cliques comprising the cut.
    std::vector<Tooth::Ptr> teeth; //<! The teeth comprising the cut.

    CliqueBank *source_bank; //<! The CliqueBank for dereferencing the cliques.
    ToothBank *source_toothbank; //<! The ToothBank for the teeth.
};

/// The external storage of a collection of HyperGraph cuts in a Relaxation.
class ExternalCuts {
public:
    /// Construct ExternalCuts with reference tour and perm.
    ExternalCuts(const std::vector<int> &tour, const std::vector<int> &perm);

    /// Add a Concorde cut.
    void add_cut(const CCtsp_lpcut_in &cc_lpcut,
                 const std::vector<int> &current_tour);

    /// Add a simple DP cut.
    void add_cut(const dominoparity &dp_cut, const double rhs,
                 const std::vector<int> &current_tour);

    /// Delete a specified set of cuts.
    void del_cuts(const std::vector<int> &delset);

    /// Return a cut corresponding to a row number index from the lp.
    const HyperGraph &get_cut(int lp_rownum) const {
        return cuts[lp_rownum - node_count];
    }

    /// Get a constant reference to the list of cuts.
    const std::vector<HyperGraph> &get_cuts() const { return cuts; }

    /// Get a constant reference to the CliqueBank.
    const CliqueBank &get_cbank() const { return clique_bank; };

    /// Get a constant reference to the ToothBank. */
    const ToothBank &get_tbank() const { return tooth_bank; };

    /// Get the column associated with an edge to be added to the lp.
    void get_col(const int end0, const int end1,
                 std::vector<int> &cmatind, std::vector<double> &cmatval);

    /// Retrieve exact/estimated duals for use in computing reduced costs.
    void get_duals(const LP::Relaxation &relax,
                   std::vector<double> &node_pi,
                   std::vector<double> &node_pi_est,
                   std::vector<double> &cut_pi,
                   std::unordered_map<Clique, double> &clique_pi) const;

private:
    /// Number of nodes in the Instance being tracked.
    /// Used to compute offsets for indices of cuts from LP::Relaxation.
    const int node_count; 

    CliqueBank clique_bank; //!< Bank for adding and dispensing cliques.
    ToothBank tooth_bank; //!< Bank for adding and dispensing teeth.

    std::vector<HyperGraph> cuts; //<! List of the cuts in the LP::Relaxation.
};

}
}

#endif
