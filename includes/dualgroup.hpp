/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief LP dual solutions.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_DUALGROUP_H
#define CMR_DUALGROUP_H

#include "lp_interface.hpp"
#include "hypergraph.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <unordered_map>

namespace CMR {
namespace LP {

/// Class template for storing dual LP solutions.
/// @tparam numtype the number type used for the vectors/pi values.
template <typename numtype>
struct DualGroup {
    DualGroup() = default; //!< Construct an empty DualGroup.

    /// Construct DualGroup from a Relaxation and HyperGraph collection.
    DualGroup(bool remove_neg, const LP::Relaxation &relax,
              const Sep::ExternalCuts &ext_cuts);

    DualGroup(DualGroup &&D) noexcept :
        node_pi(std::move(D.node_pi)),
        node_pi_est(std::move(D.node_pi_est)),
        cut_pi(std::move(D.cut_pi)),
        clique_pi(std::move(D.clique_pi)) {}

    DualGroup &operator=(DualGroup &&D) noexcept
        {
            node_pi = std::move(D.node_pi);
            node_pi_est = std::move(D.node_pi_est);
            cut_pi = std::move(D.cut_pi);
            clique_pi = std::move(D.clique_pi);
            return *this;
        }

    std::vector<numtype> node_pi; //!< Dual values for degree constraints.
    std::vector<numtype> node_pi_est; //!< Overestimates of node_pi.
    std::vector<numtype> cut_pi; //!< Dual values for cuts.

    /// Dual values/multiplicities for Cliques.
    std::unordered_map<Sep::Clique, numtype> clique_pi;
};

/**
 * Construct a DualGroup by querying the LP solver for dual values and using
 * the collection of cuts, cliques, and teeth to generate clique
 * multiplicities and node pi estimates.
 * @tparam numtype the numeric representation. Should be one of double or
 * util::Fixed64.
 * @param remove_neg if true, negative dual values will be replaced with zero.
 * @param relax the Relaxation which will be used to grab node and cut pi's.
 * @param ext_cuts the HyperGraph list of cuts in the Relaxation.
 */
template <typename numtype>
DualGroup<numtype>::DualGroup(bool remove_neg, const LP::Relaxation &relax,
                              const Sep::ExternalCuts &ext_cuts) try
{
    using HyperGraph = Sep::HyperGraph;
    using CutType = HyperGraph::Type;
    using std::vector;
    using std::cout;
    using std::cerr;
    using std::endl;

    vector<double> full_pi;
    vector<double> d_node_pi;
    vector<double> d_cut_pi;

    const vector<HyperGraph> &cuts = ext_cuts.get_cuts();
    const Sep::CliqueBank &clique_bank = ext_cuts.get_cbank();
    const vector<int> &def_tour = clique_bank.ref_tour();

    int node_count = def_tour.size();

    if (relax.num_rows() != node_count + cuts.size()) {
        cerr << "Relaxation row count: " << relax.num_rows() << ", "
             << "ExternalCuts expects: "
             << (node_count + cuts.size()) << "\n";
        throw std::logic_error("Size mismatch.");
    }

    clique_pi.reserve(clique_bank.size());
    relax.get_pi(full_pi, 0, relax.num_rows() - 1);

    node_pi = vector<numtype>(full_pi.begin(), full_pi.begin() + node_count);
    node_pi_est = node_pi;

    if (!cuts.empty())
        cut_pi = vector<numtype>(full_pi.begin() + node_count, full_pi.end());

    if (node_pi.size() != node_count || cut_pi.size() != cuts.size())
        throw std::logic_error("Node pi or cut pi size mismatch");

    if (remove_neg) {
        for (int i = 0; i < cuts.size(); ++i) {
            if (cuts[i].get_sense() == 'G') {
                if (cut_pi[i] < 0) {
                    cout << "\tCorrection: setting >= dual "
                         << cut_pi[i] << " to zero.\n";
                    cut_pi[i] = 0;
                }
            } else if (cuts[i].get_sense() == 'L') {
                if (cut_pi[i] > 0) {
                    cout << "\tCorrection: setting <= dual "
                         << cut_pi[i] << " to zero.\n";
                    cut_pi[i] = 0;
                }
            }
        }
    }

    //get clique_pi for non-domino cuts
    for (int i = 0; i < cuts.size(); ++i) {
        const HyperGraph &H = cuts[i];

        if (H.cut_type() == CutType::Non)
            throw std::logic_error("Tried to get_duals with Non cut present.");

        if (H.cut_type() == CutType::Domino)
            continue;

        numtype pival = cut_pi[i];

        for (const Sep::Clique::Ptr &clq_ref : H.get_cliques()) {
            if (clique_pi.count(*clq_ref) == 0)
                clique_pi[*clq_ref] = 0.0;

            clique_pi[*clq_ref] += pival;
        }
    }
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
        const HyperGraph &H = cuts[i];
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
} catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    throw std::runtime_error("Problem in DualGroup constructor.");
}


}
}


#endif
