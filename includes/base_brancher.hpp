/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Abstract branching control.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef CMR_BASE_BRANCHER_H
#define CMR_BASE_BRANCHER_H

#include "exec_branch.hpp"
#include "datagroups.hpp"
#include "active_tour.hpp"
#include "branch_util.hpp"

#include <list>
#include <vector>

namespace CMR {
namespace ABC {

/// Abstract base class for implementing a branching node selection rule.
class BaseBrancher {
public:
    BaseBrancher(const Data::Instance &inst, const Data::BestGroup &bestdata,
                 const Graph::CoreGraph &coregraph, LP::CoreLP &core);

    virtual ~BaseBrancher() = default;

    /// Return the next subproblem to be examined.
    virtual BranchHistory::iterator next_prob() = 0;

    /// Split the current problem, adding the subproblems to branch_history.
    void split_prob(BranchHistory::iterator &current);

    /// Branch on \p B, enforcing its constraint and instating its tour.
    void do_branch(const BranchNode &B);

    /// Unbranch on \p B and all applicable ancestors to prep next problem.
    void do_unbranch(const BranchNode &B);

    const BranchHistory &get_history() { return branch_history; }

    int verbose = 0;

protected:
    /// Adding child subproblems to branch_history.
    /// @param prob_array the pair of child subproblems to be added.
    /// This function shall be implemented to add child subproblems to the
    /// list of problems to be processed in a way that preserves the ordering
    /// criteria of the node selection rule.
    virtual void enqueue_split(BranchNode::Split prob_array) = 0;

    /// Execute variable changes if \p done was just done and \p next is next.
    void common_prep_next(const BranchNode &done, const BranchNode &next);

    /// Set next_itr to the next subproblem to be examined.
    /// This is effectively the implementation of the node selection rule.
    virtual void fetch_next() = 0;

    const Data::Instance &instance;
    const Data::BestGroup &best_data;
    const Graph::CoreGraph &core_graph;
    LP::CoreLP &core_lp;

    BranchTourFind btour_find;

    Executor exec;
    BranchHistory branch_history;

    BranchHistory::iterator next_itr;
};

}
}
#endif
