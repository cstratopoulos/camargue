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
#include <vector>

namespace CMR {
namespace ABC {

/// Depth-first search branching.
/// A child of the current node is always the next subproblem. In particular,
/// starting from the root, this class will always examine the child subproblem
/// that affirms the current best tour.
class DFSbrancher : public BaseBrancher {
public:
    DFSbrancher(const Data::Instance &inst, const LP::ActiveTour &activetour,
                const Data::BestGroup &bestdata,
                const Graph::CoreGraph &coregraph, LP::CoreLP &core);

    BranchHistory::iterator next_prob();

protected:
    void enqueue_split(BranchNode::Split prob_array);
};

}
}

#endif
