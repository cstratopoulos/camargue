#ifndef CMR_HYPERGRAPH_H
#define CMR_HYPERGRAPH_H

#include "cut_structs.hpp"
#include "lp_interface.hpp"
#include "cliq.hpp"
#include "price_util.hpp"
#include "err_util.hpp"

#include <iostream>
#include <map>
#include <stdexcept>
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

inline std::ostream &operator<<(std::ostream &os, HyperGraph::Type t)
{
    using T = HyperGraph::Type;
    if (t == T::Domino)
        os << "Domino";
    else if (t == T::Subtour)
        os << "Subtour";
    else if (t == T::Comb)
        os << "Blossom or comb";
    else
        os << "Error unkown type";
    return os;
}

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
                 std::vector<int> &cmatind, std::vector<double> &cmatval) const;

    /// Retrieve exact/estimated duals for use in computing reduced costs.
    template <typename numtype>
    void get_duals(const LP::Relaxation &relax,
                   std::vector<numtype> &node_pi,
                   std::vector<numtype> &node_pi_est,
                   std::vector<numtype> &cut_pi,
                   std::unordered_map<Clique, numtype> &clique_pi) const;

private:
    /// Number of nodes in the Instance being tracked.
    /// Used to compute offsets for indices of cuts from LP::Relaxation.
    const int node_count; 

    CliqueBank clique_bank; //!< Bank for adding and dispensing cliques.
    ToothBank tooth_bank; //!< Bank for adding and dispensing teeth.

    std::vector<HyperGraph> cuts; //<! List of the cuts in the LP::Relaxation.
};

//////////////////// TEMPLATE METHOD IMPLEMENTATION ///////////////////////////

/**
 * This method will use the LP::Relaxation to query the lp solver for dual 
 * values, and use the ExternalCuts collection of cuts, cliques, and teeth to
 * generate clique multiplicities and node pi estimates for use in exactly
 * pricing edges not currently in the Relaxation.
 * @param[in] relax the core lp relaxation, for use in grabbing node and cut
 * pi values from the lp solver.
 * @param[in,out] node_pi a vector for storing dual values associated to 
 * degree constraints
 * @param[in,out] node_pi_est a vector for storing estimated degree constraint
 * duals, computed Standard HyperGraph duals and by overestimating the 
 * Domino HyperGraph duals.
 * @param[in,out] cut_pi the dual values associated to cuts in the Relaxation.
 * @param[in,out] clique_pi a hash table of multiplicities associated to 
 * Cliques from the source CliqueBank.
 */
template <typename numtype>
void ExternalCuts::get_duals(const LP::Relaxation &relax,
                             std::vector<numtype> &node_pi,
                             std::vector<numtype> &node_pi_est,
                             std::vector<numtype> &cut_pi,
                             std::unordered_map<Clique,
                             numtype> &clique_pi) const
{
    using CutType = HyperGraph::Type;
    using std::vector;
    using std::cout;
    using std::cerr;
    using std::endl;
    
    std::runtime_error err("Problem in ExternalCuts::get_duals.");

    vector<double> d_node_pi;
    vector<double> d_cut_pi;


    if (relax.num_rows() != node_count + cuts.size()) {
        cerr << "Relaxation row count: " << relax.num_rows() << ", "
             << "ExternalCuts expects: "
             << (node_count + cuts.size()) << "\n";
        throw std::logic_error("Size mismatch in ExternalCuts::get_duals.");
    }

    clique_pi.clear();

    try {
        clique_pi.reserve(clique_bank.size());

        relax.get_pi(d_node_pi, 0, node_count - 1);

        if (!cuts.empty()) 
            relax.get_pi(d_cut_pi, node_count, relax.num_rows() - 1);
    } CMR_CATCH_PRINT_THROW("getting/allocating double pi containers.", err);
    
    try {
        node_pi = vector<numtype>(d_node_pi.begin(), d_node_pi.end());
        cut_pi = vector<numtype>(d_cut_pi.begin(), d_cut_pi.end());
        node_pi_est = node_pi;
    } CMR_CATCH_PRINT_THROW("copying to result vectors", err);

    //get clique_pi for non-domino cuts
    for (int i = 0; i < cuts.size(); ++i) {
        const Sep::HyperGraph &H = cuts[i];
        if (H.cut_type() == CutType::Domino)
            continue;

        numtype pival = cut_pi[i];

        for (const Sep::Clique::Ptr &clq_ref : H.get_cliques()) {
            if (clique_pi.count(*clq_ref) == 0)
                clique_pi[*clq_ref] = 0.0;
            
            clique_pi[*clq_ref] += pival;
        }
    }

    const vector<int> &def_tour = clique_bank.ref_tour();

    //use clique_pi to build node_pi for all cliques with nonzero pi
    for (const std::pair<Sep::Clique, numtype> &kv : clique_pi) {
        const Sep::Clique &clq = kv.first;
        numtype pival = kv.second;

        if (pival > 0.0) {
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];
                    
                    node_pi[node] += pival;
                    node_pi_est[node] += pival;
                }
            }            
        } else if (pival < 0.0) {
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];
                    
                    node_pi[node] += pival;
                }
            }
        }
    }

    //now get node_pi_est for domino cuts, skipping standard ones
    for (int i = 0; i < cuts.size(); ++i) {
        const Sep::HyperGraph &H = cuts[i];
        if (H.cut_type() != CutType::Domino)
            continue;

        numtype pival = cut_pi[i];
        if (pival <= 0.0)
            continue;
        
        const Sep::Clique::Ptr &handle_ref = H.get_cliques()[0];
        for (const Segment &seg : handle_ref->seg_list()) {
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                node_pi_est[node] += pival;
            }
        }

        for (const Sep::Tooth::Ptr &T : H.get_teeth())
            for (const Sep::Clique &tpart : T->set_pair())
                for (const Segment &seg : tpart.seg_list())
                    for (int k = seg.start; k <= seg.end; ++k) {
                        int node = def_tour[k];

                        node_pi_est[node] += pival;
                    }
    }    
}

}
}

#endif
