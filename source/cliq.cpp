#include "cliq.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::cout;
using std::cerr;

using std::unordered_map;
using std::vector;

using std::runtime_error;
using std::logic_error;
using std::exception;

using lpclique = CCtsp_lpclique;
using lpcut_in = CCtsp_lpcut_in;



namespace CMR {
namespace Sep {


/**
 * If \p cc_cliq is the clique and \p current_tour was the resident
 * best tour when the clique was found, construct a Clique where the
 * nodes in \p cc_cliq are represented as indices from \p saved_tour with
 * corresponding permutation vector \p saved_perm.
 */
Clique::Clique(const lpclique &cc_clq,
               const vector<int> &saved_tour, const vector<int> &saved_perm,
               const vector<int> &current_tour)
try
{
    for (int i = 0; i < cc_clq.segcount; ++i) {
        CMR::Segment seg(cc_clq.nodes[i].lo, cc_clq.nodes[i].hi);
        bool range_agrees = true;

        for (int k = seg.start; k <= seg.end; ++k) {
            if (saved_tour[k] != current_tour[k]) {
                range_agrees = false;
                break;
            }
        }

        if (range_agrees) {
            seglist.push_back(seg);
            continue;
        }

        vector<int> seg_nodes;
        seg_nodes.reserve(seg.size());

        for (int k = seg.start; k <= seg.end; ++k)
            seg_nodes.push_back(saved_perm[current_tour[k]]);

        std::sort(seg_nodes.begin(), seg_nodes.end());

        int k = 0;

        while (k < seg_nodes.size()) {
            int low = seg_nodes[k];

            while ((k < (seg_nodes.size() - 1)) &&
                   (seg_nodes[k + 1] == (seg_nodes[k] + 1)))
                ++k;

            seglist.push_back(CMR::Segment(low, seg_nodes[k++]));
        }
    }

    std::sort(seglist.begin(), seglist.end(), std::greater<CMR::Segment>());
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Clique CCtsp_lpclique constructor failed.");
}


/**
 * @param[in] start the start index of the Clique.
 * @param[in] end the end index of the Clique.
 * @param[in] saved_tour the indices which will be used to express the Clique.
 * @param[in] saved_perm the permutation vector for saved_tour.
 * @param[in] current_tour the tour resident when the Clique was found.
 * Constructs the Clique corresponding to the nodes
 * `current_tour[start]` up to `current_tour[end]`,
 * as indices from \p saved_tour.
 */
Clique::Clique(int start, int end,
               const vector<int> &saved_tour, const vector<int> &saved_perm,
               const vector<int> &current_tour) try
{
    bool range_agrees = true;
    CMR::Segment seg(start, end);

    for (int k = start; k <= end; ++k)
        if (saved_tour[k] != current_tour[k]) {
            range_agrees = false;
            break;
        }

    if (range_agrees) {
        seglist.push_back(seg);
        return;
    }

    vector<int> seg_nodes;
    seg_nodes.reserve(seg.size());

    for (int k = seg.start; k <= seg.end; ++k)
        seg_nodes.push_back(saved_perm[current_tour[k]]);

    std::sort(seg_nodes.begin(), seg_nodes.end());

    int k = 0;

    while (k < seg_nodes.size()) {
        int low = seg_nodes[k];

        while ((k < (seg_nodes.size() - 1)) &&
               (seg_nodes[k + 1] == (seg_nodes[k] + 1)))
            ++k;

        seglist.push_back(CMR::Segment(low, seg_nodes[k++]));
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Clique seg constructor failed.");
}


/**
 * @param[in] nodes a list of nodes in the graph as per some absolute order
 * not dependent on the current tour.
 * @param[in] perm the clique will be built using indices from perm, hence
 * implicitly reprsented in terms of the tour corresponding to \p perm.
 */
Clique::Clique(std::vector<int> &nodes, const std::vector<int> &perm)
try
{
    if (nodes.empty())
        throw logic_error("Tried to construct empty clique.");

    std::sort(nodes.begin(), nodes.end(),
              [&perm](int n1, int n2) -> bool {
                  return perm[n1] < perm[n2];
              });

    int i = 0;

    while (i < nodes.size()) {
        int low = perm[nodes[i]];

        while ((i < (nodes.size() - 1)) &&
               (perm[nodes[i + 1]] == (perm[nodes[i]] + 1)))
            ++i;

        seglist.push_back(CMR::Segment(low, perm[nodes[i++]]));
    }

    std::sort(seglist.begin(), seglist.end(), std::greater<CMR::Segment>());
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Clique nodelist constructor failed.");
}

/**
 * @param[in] saved_tour the tour that was active when this Clique was
 * constructed.
 * @returns  a vector of the literal nodes
 * obtained by dereferencing \p saved_tour for the ranges specified
 * in the segment list.
 */
vector<int> Clique::node_list(const vector<int> &saved_tour) const
{
    vector<int> result;

    for (const CMR::Segment &seg : seglist)
        for (int k = seg.start; k <= seg.end; ++k)
            result.push_back(saved_tour[k]);

    return result;
}


/**
 * @param[in] T the SimpleTooth to be represented.
 * @param[in] saved_tour the dereferencing tour for this Tooth
 * @param[in] saved_perm the dereferencing perm for this Tooth
 * @param[in] current_tour the active tour when \p T was found.
 * Uses current_tour along with saved_tour and saved_perm to convert a
 * SimpleTooth into a tooth which can be stored as a pair of Cliques.
 */
Tooth::Tooth(const SimpleTooth &T,
             const vector<int> &saved_tour, const vector<int> &saved_perm,
             const vector<int> &current_tour) try
{
  sets[0] = Clique(T.root, T.root, saved_tour, saved_perm, current_tour);
  sets[1] = Clique(T.body_start, T.body_end, saved_tour, saved_perm,
		   current_tour);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Tooth constructor failed.");
}

CliqueBank::CliqueBank(const vector<int> &tour, const vector<int> &perm)
try : saved_tour(tour), saved_perm(perm) {} catch (const exception &e) {
    throw runtime_error("CliqueBank constructor failed.");
}

Clique::Ptr CliqueBank::add_clique(const Clique &clq)
{
    if (bank.count(clq) == 0)
        bank[clq] = std::make_shared<Clique>(clq);

    return bank[clq];
}

/**
 * Same as the overload taking a Clique, but will construct a Clique
 * from the CCtsp_lpclique \p cc_cliq, with \p tour as the tour active
 * when \p cc_cliq was obtained.
 */
Clique::Ptr CliqueBank::add_clique(const CCtsp_lpclique &cc_clq,
                                    const vector<int> &tour)
{
    return add_clique(Clique(cc_clq, saved_tour, saved_perm, tour));
}

/**
 * Constructs a Clique in place and adds it to the bank, using the Clique
 * constructor taking start and end points.
 */
Clique::Ptr CliqueBank::add_clique(int start, int end, const vector<int> &tour)
{
    return add_clique(Clique(start, end, saved_tour, saved_perm, tour));
}

/**
 * @param[in] nodes shall be a list of nodes to be included in the
 * Clique.
 * @warning The elements of \p nodes are sorted by this function, but
 * unchanged otherwise.
 */
Clique::Ptr CliqueBank::add_clique(vector<int> &nodes)
{
    return add_clique(Clique(nodes, saved_perm));
}

/**
 * Steals a Clique from a bank and puts it into this CliqueBank, decrementing
 * its use count in the bank it is moved from. The reference to the Clique
 * will now be owned by this bank.
 * @param[in/out] clq_ptr the value `*clq_ptr` is inserted into this bank.
 * Then \p clq_ptr will point to `*clq_ptr` as a value in this bank, rather
 * than \p from_bank.
 * @param[in/out] from_bank the source bank that \p clq_ptr is being taken
 * from. This function will call `from_bank.del_clique(clq_ptr)`, decrementing
 * its use count in \p from_bank, possibly removing it from \p from_bank
 * entirely.
 */
void CliqueBank::steal_clique(Clique::Ptr &clq_ptr, CliqueBank &from_bank)
{
    if (!clq_ptr)
        throw logic_error("Called steal_clique on null Clique");

    Clique::Ptr new_clq = add_clique(*clq_ptr);
    from_bank.del_clique(clq_ptr);
    clq_ptr = std::move(new_clq);
}

/**
 * The Clique pointed to by \p clq_ptr will be nullified, thereby
 * decrementing the reference count of every other reference to it.
 * If its reference count in the CliqueBank drops to one, it will be
 * erased, decreasing the size of the CliqueBank.
 */
void CliqueBank::del_clique(Clique::Ptr &clq_ptr)
{
    if (!clq_ptr)
        return;

    Clique &clq = *clq_ptr;

    CliqueBank::Itr find_it = bank.find(clq);

    if (find_it == bank.end())
        return;

    clq_ptr.reset();

    if (find_it->second.use_count() <= 1)
        bank.erase(find_it);
}

ToothBank::ToothBank(const vector<int> &tour, const vector<int> &perm)
try : saved_tour(tour), saved_perm(perm) {} catch (const exception &e) {
    throw runtime_error("ToothBank constructor failed.");
}

ToothBank::ToothBank(const CliqueBank &cbank)
try : saved_tour(cbank.ref_tour()), saved_perm(cbank.ref_perm()) {}
catch (const exception &e) {
    throw runtime_error("ToothBank constructor from CliqueBank failed.");
}

/**
 * @param[in] T the SimpleTooth to add to the bank.
 * @param[in] tour the tour that was active when \p T was found.
 */
Tooth::Ptr ToothBank::add_tooth(const SimpleTooth &T, const vector<int> &tour)
{
    Tooth t(T, saved_tour, saved_perm, tour);

    if (bank.count(t) == 0)
        bank[t] = std::make_shared<Tooth>(t);

    return bank[t];
}

void ToothBank::del_tooth(Tooth::Ptr &T_ptr)
{
    if (!T_ptr)
        return;

    Tooth &T = *T_ptr;

    ToothBank::Itr find_it = bank.find(T);

    if (find_it == bank.end())
        return;

    T_ptr.reset();

    if (find_it->second.use_count() <= 1)
        bank.erase(find_it);
}

}
}
