/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Representing cuts outside the LP solver.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_HYPERGRAPH_H
#define CMR_HYPERGRAPH_H

#include "cut_structs.hpp"
#include "lp_interface.hpp"
#include "cliq.hpp"
#include "lp_util.hpp"
#include "price_util.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace CMR {
namespace Sep {

/// External representation of a cut added to the lp relaxation.
class HyperGraph {
public:
    /// Default construct an empty Type::Non HyperGraph cut.
    HyperGraph()
        : sense('\0'), source_bank(nullptr), source_toothbank(nullptr) {}

    /// Construct a HyperGraph from a Concorde cut.
    HyperGraph(CliqueBank &bank,
               const CCtsp_lpcut_in &cc_lpcut,
               const std::vector<int> &tour);

    /// Construct a HyperGraph from a simple tooth inequality.
    HyperGraph(CliqueBank &bank, ToothBank &tbank,
               const dominoparity &dp_cut, double _rhs,
               const std::vector<int> &tour);

    /// Construct a HyperGraph from ex_blossom handle/tooth indices
    HyperGraph(CliqueBank &bank,
               const std::vector<int> &blossom_handle,
               const std::vector<std::vector<int>> &tooth_edges);

    HyperGraph(HyperGraph &&H) noexcept;
    HyperGraph &operator=(HyperGraph &&H) noexcept;

    /// Change the source_bank to \p new_source_bank.
    void transfer_source(CliqueBank &new_source_bank);

    /// Turn this cut into a CCtsp_lpcut_in.
    CCtsp_lpcut_in to_lpcut_in(const std::vector<int> &active_perm,
                               bool with_skeleton) const;

    ~HyperGraph(); //!< Destruct and decrement/delete Clique/Tooth refs.

    /// Enumeration for the types of HyperGraph inequalities.
    enum Type {
        Domino = 0, //!< A domino parity inequality.
        Subtour = 1, //!< An SEC.
        Comb = 2, //!< A comb-like constraint.
        Non = 3, //!< Non HyperGraph cut: Gomory cut or branching constraint.
        Branch = 3, //!< A constraint to enforce branching.
    };

    Type cut_type() const;

    int tour_age() const { return t_age; }
    int piv_age() const { return p_age; }

    bool fresh_cut() const { return t_age <= 0 || p_age <= 0; }

    /// Get the coefficient of an edge specified by endpoints.
    double get_coeff(int end0, int end1) const;

    /// Get sparse coefficient row for a list of endpoints.
    template <typename EndPt_type>
    void get_coeffs(const std::vector<EndPt_type> &edges,
                    std::vector<int> &rmatind,
                    std::vector<double> &rmatval) const;

    char get_sense() const { return sense; }
    double get_rhs() const { return rhs; }

    const std::vector<Clique::Ptr> &get_cliques() const { return cliques; }
    const std::vector<Tooth::Ptr> &get_teeth() const { return teeth; }

    friend class ExternalCuts;

private:
    char sense; //!< The inequality sense of the cut.
    double rhs; //!< The righthand-side of the cut.

    std::vector<Clique::Ptr> cliques; //!< The cliques comprising the cut.
    std::vector<Tooth::Ptr> teeth; //!< The teeth comprising the cut.

    CliqueBank *source_bank; //!< The CliqueBank for dereferencing the cliques.
    ToothBank *source_toothbank; //!< The ToothBank for the teeth.

    int t_age;
    int p_age;
};

inline std::ostream &operator<<(std::ostream &os, HyperGraph::Type t)
{
    using T = HyperGraph::Type;
    if (t == T::Domino)
        os << "Domino";
    else if (t == T::Subtour)
        os << "Subtour";
    else if (t == T::Comb)
        os << "Blossom/Comb/Wild Thing";
    else if (t == T::Non)
        os << "Non HyperGraph Gomory cut or branch constraint.";
    return os;
}

/// The external storage of a collection of HyperGraph cuts in a Relaxation.
class ExternalCuts {
public:
    /// Construct ExternalCuts with reference tour and perm.
    ExternalCuts(const std::vector<int> &tour, const std::vector<int> &perm);

    ~ExternalCuts();

    /// Add a Concorde cut.
    void add_cut(const CCtsp_lpcut_in &cc_lpcut,
                 const std::vector<int> &current_tour);

    /// Add a simple DP cut.
    void add_cut(const dominoparity &dp_cut, const double rhs,
                 const std::vector<int> &current_tour);

    /// Add an ex_blossom cut
    void add_cut(const std::vector<int> &blossom_handle,
                 const std::vector<std::vector<int>> &tooth_edges);

    /// Add a HyperGraph cut from a pool.
    void add_cut(HyperGraph &H);

    /// Add a Non HyperGraph cut.
    void add_cut();

    void reset_ages(); //!< Reset the ages of all cuts to zero.
    void tour_age_cuts(std::vector<double> duals); //!< Update tour ages.
    void piv_age_cuts(std::vector<double> duals); //!< Update pivot ages.

    /// Delete a specified set of cuts or move them to the cut pool.
    void del_cuts(std::vector<int> &delset);

    /// Return a cut corresponding to a row number index from the lp.
    const HyperGraph &get_cut(int lp_rownum) const {
        return cuts[lp_rownum - node_count];
    }

    const std::vector<HyperGraph> &get_cuts() const { return cuts; }

    int cut_count() const { return cuts.size(); }

    int pool_count() const { return cc_pool ? cc_pool->cutcount : 0; }

    const CliqueBank &get_cbank() const { return clique_bank; }
    const ToothBank &get_tbank() const { return tooth_bank; };

    /// Get the column associated with an edge to be added to the lp.
    void get_col(int end0, int end1, std::vector<int> &cmatind,
                 std::vector<double> &cmatval) const;

    friend class Separator;

private:
    void pool_add(HyperGraph H);

    /// Number of nodes in the Instance being tracked.
    /// Used to compute offsets for indices of cuts from LP::Relaxation.
    const int node_count;

    CliqueBank clique_bank; //!< Bank for adding and dispensing cliques.
    ToothBank tooth_bank; //!< Bank for adding and dispensing teeth.

    CliqueBank pool_cliques; //!< Bank of cliques for cuts in cut_pool.

    std::vector<HyperGraph> cuts; //!< List of the cuts in the CoreLP.

    std::vector<HyperGraph> cut_pool; //!< Pool of cuts pruned from CoreLP.

    CCtsp_lpcuts *cc_pool; //!< Concorde rep of cut pool.

    CCtsp_cuttree tightcuts; //!< Cut tree for separation routines.

};

//////////////////// TEMPLATE METHOD IMPLEMENTATIONS //////////////////////////

/**
 * @tparam EndPt_type a structure derived from EndPts that stores an edge `e`
 * as a length-two array accessed as `e.end[0]` and `e.end[1]`
 * @param[in] edges the list of edges for which to generate coefficients.
 * @param[in,out] rmatind the indices of edges with nonzero coefficients.
 * @param[in,out] rmatval the coefficients corresponding to entries of
 * \p rmatind.
 */
template <typename EndPt_type>
void HyperGraph::get_coeffs(const std::vector<EndPt_type> &edges,
                            std::vector<int> &rmatind,
                            std::vector<double> &rmatval) const
{
    using std::vector;

    if (cut_type() == Type::Non)
        throw std::runtime_error("Tried HyperGraph::get_coeffs on Non cut.");

    rmatind.clear();
    rmatval.clear();

    const vector<int> &def_tour = source_bank->ref_tour();
    int ncount = def_tour.size();

    std::map<int, int> coeff_map;
    vector<bool> node_marks(ncount, false);

    if (cut_type() != Type::Domino) {
        for (const Clique::Ptr &clq_ref : cliques) {
            for (const Segment &seg : clq_ref->seg_list())
                for (int k = seg.start; k <= seg.end; ++k)
                    node_marks[def_tour[k]] = true;

            for (int i = 0; i < edges.size(); ++i) {
                if (node_marks[edges[i].end[0]] !=
                    node_marks[edges[i].end[1]]) {
                    if (coeff_map.count(i))
                        coeff_map[i] += 1;
                    else
                        coeff_map[i] = 1;
                }
            }

            node_marks = vector<bool>(ncount, false);
        }

        rmatind.reserve(coeff_map.size());
        rmatval.reserve(coeff_map.size());

        for (std::pair<const int, int> &kv : coeff_map) {
            rmatind.push_back(kv.first);
            rmatval.push_back(kv.second);
        }

        return;
    } //else it is a domino cut

    const Clique::Ptr &handle_ref = cliques[0];
    for (const Segment &seg : handle_ref->seg_list())
        for (int k = seg.start; k <= seg.end; ++k)
            node_marks[def_tour[k]] = true;

    for (int i = 0; i < edges.size(); ++i) {
        bool e0 = node_marks[edges[i].end[0]];
        bool e1 = node_marks[edges[i].end[1]];

        if (e0 && e1) {
            if (coeff_map.count(i))
                coeff_map[i] += 2;
            else
                coeff_map[i] = 2;
        } else if (e0 != e1) {
            if (coeff_map.count(i))
                coeff_map[i] += 1;
            else
                coeff_map[i] = 1;
        }
    }

    node_marks = vector<bool>(ncount, false);

    const vector<int> &tooth_perm = source_toothbank->ref_perm();

    for (const Tooth::Ptr &t_ref : teeth) {
        const Clique &root_clq = t_ref->set_pair()[0];
        const Clique &bod_clq = t_ref->set_pair()[1];

        for (const Segment &seg : bod_clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k)
                node_marks[def_tour[k]] = true;

        for (int i = 0; i < edges.size(); ++i) {
            int e0 = edges[i].end[0];
            int e1 = edges[i].end[1];

            bool bod_e0 = node_marks[e0];
            bool bod_e1 = node_marks[e1];

            if (bod_e0 && bod_e1) {
                if (coeff_map.count(i))
                    coeff_map[i] += 2;
                else
                    coeff_map[i] = 2;
            } else if (bod_e0 != bod_e1) {
                if (bod_e0) {
                    if (root_clq.contains(tooth_perm[e1])) {
                        if (coeff_map.count(i))
                            coeff_map[i] += 1;
                        else coeff_map[i] = 1;
                    }
                } else { //bod_e1
                    if (root_clq.contains(tooth_perm[e0])) {
                        if (coeff_map.count(i))
                            coeff_map[i] += 1;
                        else coeff_map[i] = 1;
                    }
                }
            }
        }

        std::fill(node_marks.begin(), node_marks.end(), false);
    }

    int coeff_size = 0;

    for (std::pair<const int, int> &kv : coeff_map)
        if (kv.second % 2 == 0)
            ++coeff_size;

    rmatind.reserve(coeff_size);
    rmatval.reserve(coeff_size);

    for (std::pair<const int, int> &kv : coeff_map) {
        int ind = kv.first;
        double coeff = kv.second;
        coeff /= 2;
        coeff = floor(coeff);
        if (fabs(coeff) >= Epsilon::Zero) {
            rmatind.push_back(ind);
            rmatval.push_back(coeff);
        }
    }
}

}
}

#endif
