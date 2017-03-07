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

    SplitIter next_level();

    void do_branch(BranchNode &B);
    void do_unbranch(BranchNode &B);

    std::vector<int> get_tour(const BranchNode &B);

    const BranchHistory &get_history() { return branch_history; }

private:
    BranchHistory::iterator next_prob();

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
