/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Branching execution.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_EXEC_BRANCH_H
#define CMR_EXEC_BRANCH_H

#include "active_tour.hpp"
#include "core_lp.hpp"
#include "datagroups.hpp"
#include "cliq.hpp"
#include "branch_node.hpp"
#include "lp_util.hpp"
#include "branch_util.hpp"

#include <array>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace CMR {
namespace ABC {

class Executor {
public:
    /// Construct an Executor using data from an existing solution process.
    Executor(const Data::Instance &inst, const LP::ActiveTour &activetour,
             const Data::BestGroup &bestdata,
             const Graph::CoreGraph &coregraph, LP::CoreLP &core);

    ScoreTuple branch_edge(); //!< Get the next edge to branch on.

    /// Create the children nodes of \p parent for branching on \p branch_edge.
    BranchNode::Split split_problem(ScoreTuple &branch_tuple,
                                    BranchNode &parent);

    /// Get a tour satisfying the branching in \p edge_stats.
    void branch_tour(const std::vector<EndsDir> &edge_stats,
                     const std::vector<int> &start_tour_nodes,
                     bool &found_tour, bool &feas,
                     std::vector<int> &tour, double &tour_val,
                     bool for_use);

    /// Clamp a variable as indicated by \p current_node.
    void clamp(const BranchNode &current_node);

    /// Undo the clamp done on \p current_node.
    void unclamp(const BranchNode &current_node);

private:
    const Data::Instance &instance;
    const LP::ActiveTour &active_tour;
    const Data::BestGroup &best_data;
    const Graph::CoreGraph &core_graph;

    LP::CoreLP &core_lp;

    /// Used by branch_tour to track the edges in \p edge_stats.
    /// An instance.node_count() dimensional upper triangular matrix with
    /// `inds_table(i, j)` being 'A' if the edge is to be avoided, 'W' if
    /// the edge is wanted (fixed to one), and '\0' otherwise.
    util::SquareUT<char> inds_table;

    /// Used by branch_tour to track fixed edges.
    /// An `instance.node_count()` length array initialized to all zero, and
    /// incremented for the ends of every BranchNode::Dir::Up edge in
    /// edge_stats. No tour can exist if any entry is greater than two.
    std::vector<int> fix_degrees;

    /// Used by branch_tour to track avoided edges.
    /// Like fix_degrees, but for track BranchNode::Dir::Down edges in
    /// edge_stats. Initialized to all `instance.node_count() - 1`, and
    /// decremented. No tour can exist if an entry is less than two.
    std::vector<int> avail_degrees;
};

}
}

#endif
