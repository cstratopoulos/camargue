/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Abstract branching control.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
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

class BaseBrancher {
public:
    BaseBrancher(const Data::Instance &inst, const LP::ActiveTour &activetour,
                 const Data::BestGroup &bestdata,
                 const Graph::CoreGraph &coregraph, LP::CoreLP &core);

    /// Implementation of a node selection rule.
    /// Node selection rules such as best bound, best estimate, depth first,
    /// etc shall be implemented by this function which picks the next
    /// unvisited element from the branch_history list.
    virtual BranchHistory::iterator next_prob() = 0;

    /// Split the current problem, adding the subproblems to branch_history.
    void split_prob(BranchHistory::iterator &current);

    /// Branch on \p B, enforcing its constraint and instating its tour.
    void do_branch(const BranchNode &B);

    /// Unbranch on \p B and all applicable ancestors to prep next problem.
    void do_unbranch(const BranchNode &B);

    const BranchHistory &get_history() { return branch_history; }

protected:
    /// Implementation of adding child subproblems to branch_history.
    /// This function shall be implemented to add child subproblems to the
    /// list of problems to be processed in a way that preserves the ordering
    /// criteria of the node selection rule.
    virtual void enqueue_split(BranchNode::Split prob_array) = 0;

    /// Execute variable changes if \p done was just done and \p next is next.
    void common_prep_next(const BranchNode &done, const BranchNode &next);

    const Data::Instance &instance;
    const Data::BestGroup &best_data;
    const Graph::CoreGraph &core_graph;
    LP::CoreLP &core_lp;

    Executor exec;
    BranchHistory branch_history;
};

}
}
#endif
