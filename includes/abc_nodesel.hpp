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
/// The child of the current node that agrees with the current tour is always
/// examined first.
class DFSbrancher : public BaseBrancher {
public:
    DFSbrancher(const Data::Instance &inst, const LP::ActiveTour &active_tour,
                const Data::BestGroup &best_data,
                const Graph::CoreGraph &core_graph, LP::CoreLP &core_lp);

    BranchHistory::iterator next_prob();

protected:
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
    void enqueue_split(BranchNode::Split prob_array);
};

}
}

#endif
