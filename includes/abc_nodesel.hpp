/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Classes for node selection rules.
 * This file contains class definitions for all implemented ABC node selection
 * rules. Implementations are in their own cpp files, indicated in the class
 * documentation.
 * @see CMR::ABC::BaseBrancher for the meaning of inherited virtual methods.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef CMR_ABC_NODESEL_H
#define CMR_ABC_NODESEL_H

#include "datagroups.hpp"
#include "active_tour.hpp"
#include "branch_node.hpp"
#include "base_brancher.hpp"
#include "qpref_brancher.hpp"

#include <vector>

namespace CMR {
namespace ABC {

/// Depth-first search branching.
/// The child of the current node that agrees with the current tour is always
/// examined first.
class DFSbrancher : public BaseBrancher {
public:
    DFSbrancher(const Data::Instance &inst, const Data::BestGroup &best_data,
                const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    void fetch_next();

    void enqueue_split(BranchNode::Split prob_array);
};

/// Alias declaration for "best tour" branching.
/// @see BranchNode::tour_worse for implementation and tie-breaking rules.
using TourBrancher = QprefBrancher<BranchNode::tour_worse>;

/// Alias declaration for best-first search branching.
/// @see BranchNode::bound_worse for implementation and tie-breaking rules.
using BoundBrancher = QprefBrancher<BranchNode::bound_worse>;

/// Interleaved best-estimate and best-first search branching.
/// Tour length is used as the primary node selection criterion, with a
/// best bound node selected every InterBrancher::BestFreq nodes.
class InterBrancher : public BaseBrancher {
public:
    InterBrancher(const Data::Instance &inst, const Data::BestGroup &best_data,
                  const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    void fetch_next();
    void enqueue_split(BranchNode::Split prob_array);

private:
    static constexpr int BestFreq = 10;
    int node_num = 1;

    /**@name Priority queue adaptors.
     * This class effectively uses a priority queue of problems, but the
     * actual std::priority_queue cannot be used because occasionally we
     * need to extract the lowest bound element. These functions adapt
     * std::push_heap, std::make_heap, and std::pop_heap to modify
     * prob_q on best tour calls.
     */
    ///@{

    static void heap_push(std::vector<BranchHistory::iterator> &target_q,
                          BranchHistory::iterator itr);

    static void heap_make(std::vector<BranchHistory::iterator> &target_q);

    static void heap_pop(std::vector<BranchHistory::iterator> &target_q);

    /// returns true iff \p A has a better (i.e.,  _lower_ )estimate than \p B.
    /// @remark this is the _opposite_ of BranchNode::bound_worse.
    static bool better_bound(const BranchHistory::iterator &A,
                           const BranchHistory::iterator &B)
        { return A->estimate < B->estimate; }

    ///@}

    std::vector<BranchHistory::iterator> prob_q;
};

}
}

#endif
