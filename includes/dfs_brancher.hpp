/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief DFS branching control.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef CMR_DFS_BRANCHER_H
#define CMR_DFS_BRANCHER_H

#include "exec_branch.hpp"
#include "datagroups.hpp"
#include "active_tour.hpp"
#include "branch_util.hpp"

#include <list>
#include <vector>

namespace CMR {
namespace ABC {

class DFSbrancher {
public:
    DFSbrancher(const Data::Instance &inst, const LP::ActiveTour &activetour,
                const Data::BestGroup &bestdata,
                const Graph::CoreGraph &coregraph, LP::CoreLP &core);

    BranchHistory::iterator next_prob();

    /// Split the current problem, adding the subproblems to branch_history.
    void split_prob(BranchHistory::iterator &current);

    /// Branch on \p B, enforcing its constraint and instating its tour.
    void do_branch(const BranchNode &B);

    /// Unbranch on \p B and all applicable ancestors to prep next problem.
    void do_unbranch(const BranchNode &B);

    /// Expand the branching tour for node \p B into a node list.
    std::vector<int> get_tour(const BranchNode &B);

    const BranchHistory &get_history() { return branch_history; }

private:
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
