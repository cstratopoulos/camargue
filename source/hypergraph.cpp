#include "hypergraph.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/cuttree.h>
}

using std::unordered_map;
using std::vector;
using std::pair;

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

using lpcut_in = CCtsp_lpcut_in;
using lpclique = CCtsp_lpclique;

namespace CMR {

namespace Eps = Epsilon;

namespace Sep {


/**
 * The cut corresponding to \p cc_lpcut will be represented using
 * Clique pointers from \p bank, assuming that the cut was found with
 * \p tour as the resident best tour.
 */
HyperGraph::HyperGraph(CliqueBank &bank, const lpcut_in &cc_lpcut,
                       const vector<int> &tour) try :
    sense(cc_lpcut.sense), rhs(cc_lpcut.rhs), source_bank(&bank),
    source_toothbank(nullptr),
    t_age(LP::CutAge::Babby), p_age(LP::CutAge::Babby)
{
    for (int i = 0; i < cc_lpcut.cliquecount; ++i) {
        lpclique &cc_clq = cc_lpcut.cliques[i];
        cliques.push_back(source_bank->add_clique(cc_clq, tour));
    }

} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("HyperGraph CC lpcut_in constructor failed.");
}

/**
 * The cut corresponding to \p dp_cut will be represented using Clique
 * pointers from \p bank and Tooth pointers from \p tbank, assuming the cut
 * was found with \p tour as the resident best tour. The righthand side of the
 * cut stored shall be \p _rhs.
 */
HyperGraph::HyperGraph(CliqueBank &bank, ToothBank &tbank,
                       const dominoparity &dp_cut, const double _rhs,
                       const std::vector<int> &tour) try :
    sense('L'), rhs(_rhs), source_bank(&bank), source_toothbank(&tbank),
    t_age(LP::CutAge::Babby), p_age(LP::CutAge::Babby)
{
    vector<int> nodes(dp_cut.degree_nodes);
    for(int &n : nodes)
        n = tour[n];

    cliques.push_back(source_bank->add_clique(nodes));

    for (const SimpleTooth &T : dp_cut.used_teeth)
        teeth.push_back(source_toothbank->add_tooth(T, tour));

} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("HyperGraph dominoparity constructor failed.");
}

/**
 * Constructs a HyperGraph from the handle nodes and edge teeth nodes of an
 * exact primal blossom inequality.
 * @param[in] bank the source CliqueBank.
 * @param[in] blossom_handle the handle of the blossom inequality;
 * @param[in] tooth_edges a vector of vectors of length two, consisting of
 * the end points of edges which are the blossom teeth.
 */
HyperGraph::HyperGraph(CliqueBank &bank,
                       const vector<int> &blossom_handle,
                       const vector<vector<int>> &tooth_edges) try
    : sense('G'), rhs ((3 * tooth_edges.size()) + 1), source_bank(&bank),
      source_toothbank(nullptr),
      t_age(LP::CutAge::Babby), p_age(LP::CutAge::Babby)
{
    vector<int> handle(blossom_handle);
    cliques.push_back(source_bank->add_clique(handle));

    for (vector<int> tooth_edge : tooth_edges)
        cliques.push_back(source_bank->add_clique(tooth_edge));
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("HyperGraph ex_blossom constructor failed.");
}

/**
 * The moved-from Hypergraph \p H is left in a null but valid state, as if it
 * had been default constructed.
 */
HyperGraph::HyperGraph(HyperGraph &&H) noexcept
    : sense(H.sense), rhs(H.rhs),
      cliques(std::move(H.cliques)), teeth(std::move(H.teeth)),
      source_bank(H.source_bank), source_toothbank(H.source_toothbank),
      t_age(LP::CutAge::Babby), p_age(LP::CutAge::Babby)
{
    H.sense = '\0';
    H.source_bank = nullptr;
    H.source_toothbank = nullptr;
    H.t_age = LP::CutAge::Babby;
    H.p_age = LP::CutAge::Babby;
}

/**
 * The move-assigned-from HyperGraph \p H is left null but valid as if default
 * constructed.
 */
HyperGraph &HyperGraph::operator=(Sep::HyperGraph &&H) noexcept
{
    sense = H.sense;
    H.sense = '\0';

    rhs = H.rhs;

    cliques = std::move(H.cliques);
    teeth = std::move(H.teeth);

    source_bank = H.source_bank;
    H.source_bank = nullptr;

    source_toothbank = H.source_toothbank;
    H.source_toothbank = nullptr;

    t_age = H.t_age;
    p_age = H.p_age;

    H.t_age = LP::CutAge::Babby;
    H.p_age = LP::CutAge::Babby;

    return *this;
}

/**
 * Ownership of the cliques in this HyperGraph is transferred to
 * \p new_source_bank, and \p new_source_bank becomes source_bank.
 */
void HyperGraph::transfer_source(CliqueBank &new_source_bank) try
{
    if (cut_type() == Type::Non || cut_type() == Type::Branch)
        throw logic_error("Tried to transfer_source on Non cut");

    for (Clique::Ptr &clq_ptr : cliques)
        new_source_bank.steal_clique(clq_ptr, *source_bank);

    source_bank = &new_source_bank;
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("HyperGraph::transfer_source failed.");
}

/**
 * @param active_perm the permutation vector for the active best tour.
 * @param with_skeleton true iff a skeleton should be constructed for the cut
 * as well.
 * @returns a CCtsp_lpcut_in representation of this HyperGraph cut.
 */
CCtsp_lpcut_in HyperGraph::to_lpcut_in(const vector<int> &active_perm,
                                       bool with_skeleton) const
{
    if (source_toothbank != nullptr)
        throw runtime_error("Called HyperGraph::to_lpcut_in on domino");
    if (source_bank == nullptr)
        throw runtime_error("Called HyperGraph::to_lpcut_in w no source_bank");

    runtime_error err("Problem in HyperGraph::to_lpcut_in");
    int rval = 0;
    CCtsp_lpcut_in result;

    CCtsp_init_lpcut_in(&result);
    result.rhs = rhs;
    result.sense = sense;

    int num_cliques = cliques.size();
    auto cleanup = util::make_guard([&rval, &result]
                                    {
                                        if (rval) CCtsp_free_lpcut_in(&result);
                                    });

    if ((rval = CCtsp_create_lpcliques(&result, num_cliques)))
        throw err;

    const vector<int> &def_tour = source_bank->ref_tour();
    for (int i = 0; i < num_cliques; ++i) {
        const Clique::Ptr &clq_ref = cliques[i];
        vector<int> clq_nodes;

        try { clq_nodes = clq_ref->node_list(def_tour); }
        catch (const exception &e) {
            rval = 1;
            cerr << e.what() << " getting clique nodes" << endl;
            throw err;
        }

        for (int &n : clq_nodes)
            n = active_perm[n];

        CCtsp_lpclique *cur_clq = (result.cliques + i);

        if ((rval = CCtsp_array_to_lpclique(&clq_nodes[0], clq_nodes.size(),
                                            cur_clq)))
            throw err;
    }

    if (with_skeleton) {
        if ((rval = CCtsp_construct_skeleton(&result, active_perm.size())))
            throw err;
    }

    return result;
}

HyperGraph::~HyperGraph()
{
    if (source_bank != nullptr)
        for (Clique::Ptr &ref : cliques)
            source_bank->del_clique(ref);

    if (source_toothbank != nullptr)
        for (Tooth::Ptr &ref : teeth)
            source_toothbank->del_tooth(ref);
}

HyperGraph::Type HyperGraph::cut_type() const
{
    using Type = HyperGraph::Type;

    if (source_bank == nullptr && source_toothbank == nullptr)
        return Type::Non;

    if (!teeth.empty())
        return Type::Domino;

    return (cliques.size() == 1) ? Type::Subtour : Type::Comb;
}

double HyperGraph::get_coeff(int end0, int end1) const
{
    if (end0 == end1)
        throw logic_error("Edge has same endpoints in HyperGraph::get_coeff.");

    if (cut_type() == Type::Non)
        throw logic_error("Tried HyperGraph::get_coeff on Non cut.");

    double result = 0.0;

    if (cut_type() != Type::Domino) {
        const vector<int> &perm = source_bank->ref_perm();

        int end0_ind = perm[end0];
        int end1_ind = perm[end1];

        for (const Clique::Ptr &clq_ref : cliques) {
            bool contains_end0 = clq_ref->contains(end0_ind);
            bool contains_end1 = clq_ref->contains(end1_ind);

            if (contains_end0 != contains_end1)
                result += 1.0;
        }

        return result;
    }

    //else it is a Simple DP
    int pre_result = 0;
    //get handle coeffs
    const vector<int> &handle_perm = source_bank->ref_perm();
    int end0_ind = handle_perm[end0];
    int end1_ind = handle_perm[end1];

    const Clique::Ptr &handle_clq = cliques[0];

    bool contains_end0 = handle_clq->contains(end0_ind);
    bool contains_end1 = handle_clq->contains(end1_ind);

    if (contains_end0 && contains_end1) //in E(H)
        pre_result += 2;
    else if (contains_end0 != contains_end1) // in delta(H)
        pre_result += 1;

    const vector<int> &tooth_perm = source_toothbank->ref_perm();

    end0_ind = tooth_perm[end0];
    end1_ind = tooth_perm[end1];
    contains_end0 = false;
    contains_end1 = false;

    for (const Tooth::Ptr &tooth : teeth) {
        const Clique &root_clq = tooth->set_pair()[0];
        const Clique &bod_clq = tooth->set_pair()[1];

        bool root_end0 = false;
        bool root_end1 = false;

        if ((root_end0 = root_clq.contains(end0_ind)))
            if (bod_clq.contains(end1_ind)) {
                pre_result += 1;
                continue;
            }

        if ((root_end1 = root_clq.contains(end1_ind)))
            if (bod_clq.contains(end0_ind)) {
                pre_result += 1;
                continue;
            }

        if (root_end0 || root_end1)
            continue;

        if (bod_clq.contains(end0_ind) && bod_clq.contains(end1_ind))
            pre_result += 2;
    }

    pre_result /= 2;
    return static_cast<double>(pre_result);
}

ExternalCuts::ExternalCuts(const vector<int> &tour, const vector<int> &perm)
try : node_count(tour.size()), clique_bank(tour, perm), tooth_bank(tour, perm),
      pool_cliques(tour, perm)
{
    int ncount = node_count;
    if (CCtsp_init_cutpool(&ncount, NULL, &cc_pool))
        throw runtime_error("CCtsp_init_cutpool failed");

    CCpq_cuttree_init(&tightcuts);

    if (CCpq_cuttree_trivial(&tightcuts, ncount, 0)) {
        CCtsp_free_cutpool(&cc_pool);
        throw runtime_error("CCpq_cuttree_trivial failed");
    }
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("ExternalCuts constructor failed.");
}

ExternalCuts::~ExternalCuts()
{
    CCtsp_free_cutpool(&cc_pool);
    CCpq_cuttree_freetree(&tightcuts);
}

/**
 * @param[in] cc_lpcut the Concorde cut to be added.
 * @param[in] current_tour the tour active when cc_lpcut was found.
 */
void ExternalCuts::add_cut(const lpcut_in &cc_lpcut,
                           const vector<int> &current_tour)
{
    cuts.emplace_back(clique_bank, cc_lpcut, current_tour);
}

/**
 * @param[in] dp_cut the domino parity inequality to be added.
 * @param[in] rhs the righthand-side of the cut.
 * @param[in] current_tour the tour active when dp_cut was found.
 */
void ExternalCuts::add_cut(const dominoparity &dp_cut, const double rhs,
                           const vector<int> &current_tour)
{
    cuts.emplace_back(clique_bank, tooth_bank, dp_cut, rhs, current_tour);
}

/**
 * @param[in] blossom_handle the handle of the blossom inequality;
 * @param[in] tooth_edges a vector of vectors of length two, consisting of
 * the end points of edges which are the blossom teeth.
 */
void ExternalCuts::add_cut(const vector<int> &blossom_handle,
                           const vector<vector<int>> &tooth_edges)
{
    cuts.emplace_back(clique_bank, blossom_handle, tooth_edges);
}

/**
 * This function is meant to be invoked after a call to PoolSep which takes
 * \p H from cut_pool and identifies that it is violated by the current LP
 * solution, hence adding it back to the LP and the vector cuts.
 */
void ExternalCuts::add_cut(HyperGraph &H)
{
    H.transfer_source(clique_bank);
    cuts.emplace_back(std::move(H));
}

/**
 * Add a branching constraint or Non HyperGraph cut to the list. Maintains
 * indexing that agrees with the Relaxation for bookkeeping and cut pruning
 * purposes.
 */
void ExternalCuts::add_cut() { cuts.emplace_back(); }

void ExternalCuts::reset_ages()
{
    for (HyperGraph &H : cuts) {
        H.t_age = LP::CutAge::Babby;
        H.p_age = LP::CutAge::Babby;
    }
}

void ExternalCuts::tour_age_cuts(vector<double> duals)
{
    if (duals.size() != cuts.size()) {
        cerr << "Duals size " << duals.size() << " vs " << cuts.size()
             << " cuts" << endl;
        throw runtime_error("Size mismatch in ExternalCuts::tour_age_cuts");
    }

    for (int i = 0; i < duals.size(); ++i)
        if (duals[i] >= Eps::DualDust)
            cuts[i].t_age = LP::CutAge::Babby;
        else
            cuts[i].t_age += 1;
}

void ExternalCuts::piv_age_cuts(vector<double> duals)
{
    if (duals.size() != cuts.size()) {
        cerr << "Duals size " << duals.size() << " vs " << cuts.size()
             << " cuts" << endl;
        throw runtime_error("Size mismatch in ExternalCuts::piv_age_cuts");
    }

    for (int i = 0; i < duals.size(); ++i)
        if (duals[i] >= Eps::DualDust)
            cuts[i].p_age = LP::CutAge::Babby;
        else
            cuts[i].p_age += 1;
}

/**
 * @param delset vector with nonzero entries indicating cuts to be deleted
 * from an LP::Relaxation.
 * @pre `delset[i] == 1` if `get_cut(i)` is to be deleted immediately, and
 * `delset[i] == 3` if `get_cut(i)` should be considered for addition to the
 * cut pool.
 * @post \p delset will be binary valued with all 3's changed to 1's.
 */
void ExternalCuts::del_cuts(vector<int> &delset)
{
    using CutType = HyperGraph::Type;

    for (int i = 0; i < cuts.size(); ++i) {
        HyperGraph &H = cuts[i];
        CutType Htype = H.cut_type();

        if (delset[i + node_count] == 0)
            continue;

        if (delset[i + node_count] == 3) {
            if (Htype == CutType::Domino) { //unimplemented for now
                // H.t_age = LP::CutAge::Babby;
                // H.p_age = LP::CutAge::Babby;
                // H.transfer_source(pool_cliques);
                // cut_pool.emplace_back(std::move(H));
            } else if (Htype == CutType::Comb) {
                pool_add(std::move(H));
            }
            delset[i + node_count] = 1;
        }
        H.sense = 'X';
    }

    util::erase_remove(cuts, [](const HyperGraph &H)
                       { return H.sense == 'X'; });
}


/**
 * @param[in] end0 one end of the edge to be added
 * @param[in] end1 the other end of the edge to be added
 * @param[in,out] cmatind the indices of the rows having nonzero
 * coefficients for the new edge
 * @param[in,out] cmatval the coefficients corresponding to entries of
 * \p cmatind.
 */
void ExternalCuts::get_col(const int end0, const int end1,
                           vector<int> &cmatind, vector<double> &cmatval) const
{
    runtime_error err("Problem in ExternalCuts::get_col");

    if (end0 == end1) {
        cerr << "Edge has same endpoints.\n";
        throw err;
    }

    cmatind.clear();
    cmatval.clear();

    int lp_size = node_count + cuts.size();

    try {
        cmatind.reserve(lp_size);
        cmatval.reserve(lp_size);

        cmatind.push_back(end0);
        cmatind.push_back(end1);

        cmatval.push_back(1.0);
        cmatval.push_back(1.0);

        for (int i = 0; i < cuts.size(); ++i) {
            int index = i + node_count;
            double coeff = cuts[i].get_coeff(end0, end1);

            if (coeff != 0.0) {
                cmatind.push_back(index);
                cmatval.push_back(coeff);
            }
        }
    } CMR_CATCH_PRINT_THROW("Couldn't push back column coeffs/inds", err);
}

void ExternalCuts::pool_add(HyperGraph H)
{
    runtime_error err("Problem in ExternalCuts::pool_add");

    if (H.cut_type() != HyperGraph::Type::Comb) {
        cerr << "Tried to add cut of type " << H.cut_type() << " to CC pool"
             << endl;
        throw err;
    }

    CCtsp_lpcut_in c;
    auto c_guard = util::make_guard([&c] { CCtsp_free_lpcut_in(&c); });

    try { c = H.to_lpcut_in(H.source_bank->ref_perm(), true); }
    CMR_CATCH_PRINT_THROW("getting lpcut_in from HyperGraph", err);

    if (CCtsp_add_to_cutpool_lpcut_in(cc_pool, &c))
        throw runtime_error("CCtsp_add_to_cutpool_lpcut_in failed");
}


}
}
