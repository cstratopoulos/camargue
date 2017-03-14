/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Classes for node selection rules.
 * This file contains class definitions for all implemented ABC node selection
 * rules. Implementations are in their own cpp files, indicated in the class
 * documentation.
 * @see CMR::ABC::BaseBrancher for the meaning of inherited virtual methods.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef CMR_DFS_BRANCHER_H
#define CMR_DFS_BRANCHER_H

#include "datagroups.hpp"
#include "active_tour.hpp"
#include "branch_node.hpp"
#include "base_brancher.hpp"

#include <list>
#include <queue>
#include <vector>

namespace CMR {
namespace ABC {

/// Depth-first search branching.
/// The child of the current node that agrees with the current tour is always
/// examined first.
class DFSbrancher : public BaseBrancher {
public:
    DFSbrancher(const Data::Instance &inst, const LP::ActiveTour &active_tour,
                const Data::BestGroup &best_data,
                const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    BranchHistory::iterator peek_next() { return next_prob(); }

    void enqueue_split(BranchNode::Split prob_array);
};

/// Best estimate search branching.
/// The node with best estimated tour length is examined first, as determined
/// by the tour computed by ExecBranch::branch_tour.
class TourBrancher : public BaseBrancher {
public:
    TourBrancher(const Data::Instance &inst, const LP::ActiveTour &active_tour,
                 const Data::BestGroup &best_data,
                 const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    BranchHistory::iterator peek_next();
    void enqueue_split(BranchNode::Split prob_array);

private:
    std::priority_queue<BranchHistory::iterator,
                        std::vector<BranchHistory::iterator>,
                        decltype(&BranchNode::tour_worse)> prob_q;
};

/// Best first (best bound) search branching.
/// The node with the lowest strong branch objective value estimate is
/// examined first.
class BoundBrancher : public BaseBrancher {
public:
    BoundBrancher(const Data::Instance &inst,
                  const LP::ActiveTour &active_tour,
                  const Data::BestGroup &best_data,
                  const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    BranchHistory::iterator peek_next();
    void enqueue_split(BranchNode::Split prob_array);

private:
    std::priority_queue<BranchHistory::iterator,
                        std::vector<BranchHistory::iterator>,
                        decltype(&BranchNode::bound_worse)> prob_q;
};

/// Interleaved best-estimate and best-first search branching.
/// Tour length is used as the primary node selection criterion, with a
/// best bound node selected every InterBrancher::BestFreq nodes.
class InterBrancher : public BaseBrancher {
public:
    InterBrancher(const Data::Instance &inst,
                  const LP::ActiveTour &active_tour,
                  const Data::BestGroup &best_data,
                  const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
    BranchHistory::iterator peek_next();
    void enqueue_split(BranchNode::Split prob_array);

private:
    static constexpr int BestFreq = 10;
    int node_num = 1;

    static void heap_push(std::vector<BranchHistory::iterator> &target_q,
                          BranchHistory::iterator itr);

    static void heap_make(std::vector<BranchHistory::iterator> &target_q);

    static void heap_pop(std::vector<BranchHistory::iterator> &target_q);

    /// returns true iff \p A has a better (i.e.,  _lower_ )estimate than \p B.
    /// @remark this is the _opposite_ of BranchNode::bound_worse.
    static bool better_bound(const BranchHistory::iterator &A,
                           const BranchHistory::iterator &B)
        { return A->estimate < B->estimate; }


    std::vector<BranchHistory::iterator> prob_q;
};

}
}

#endif
