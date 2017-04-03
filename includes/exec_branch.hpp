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
#include "branch_tour.hpp"
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
    Executor(const Data::Instance &inst, const Data::BestGroup &bestdata,
             const Graph::CoreGraph &coregraph, LP::CoreLP &core,
             BranchTourFind &btourfind);

    ScoreTuple branch_edge(); //!< Get the next edge to branch on.

    /// Create the children nodes of \p parent for branching on \p branch_edge.
    BranchNode::Split split_problem(ScoreTuple &branch_tuple,
                                    BranchNode &parent);

    /// Clamp a variable as indicated by \p current_node.
    void clamp(const BranchNode &current_node);

    /// Undo the clamp done on \p current_node.
    void unclamp(const BranchNode &current_node);

    int verbose = 0;

private:
    const Data::Instance &instance;
    const LP::ActiveTour &active_tour;
    const Data::BestGroup &best_data;
    const Graph::CoreGraph &core_graph;

    LP::CoreLP &core_lp;

    BranchTourFind &btour_find;

};

}
}

#endif
