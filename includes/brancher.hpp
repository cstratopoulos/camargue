/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**
 * @file
 * @brief Class for managing an ABC search.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCHER_H
#define CMR_BRANCHER_H

#include "datagroups.hpp"
#include "core_lp.hpp"
#include "graph.hpp"
#include "branch_util.hpp"

#include <array>
#include <iostream>
#include <functional>
#include <utility>
#include <vector>

namespace CMR {
namespace ABC {

/// Generation and dispensing of branching subproblems.
class Brancher {
public:
    Brancher(LP::Relaxation &lp_rel, const Graph::CoreGraph &coregraph,
             const LP::TourBasis &tbase, const double &tourlen,
             ContraStrat strat);

    std::array<Problem, 2> next_level(); //!< Go a level deeper in search tree.

    ScoreTuple next_branch_obj(); //!< Describes the next edge to branch on.

    /// Compute a tour which satisfies a branching problem.
    std::vector<int> branch_tour(const Data::Instance &inst,
                                 std::vector<int> &start_tour_nodes);

    void do_branch(Problem &prob); //!< Enforce a branching problem.
    void undo_branch(Problem &prob); //!< Unenforce a branching problem.

    using EdgeStats = std::pair<Graph::Edge, int>;

    const std::vector<EdgeStats> &branch_statuses() const
        { return branch_stats; }

private:
    LP::Relaxation &lp_relax;
    const Graph::CoreGraph &core_graph;
    const LP::TourBasis &tour_base;
    const double &tour_len;

    const ContraStrat contra_strategy;
    const std::function<void(LP::Relaxation&, int, double)> contra_enforce;
    const std::function<void(LP::Relaxation&, int, double)> contra_undo;

    std::vector<EdgeStats> branch_stats;
};

/// Enforce a Contra branch on an edge by fixing bounds.
void contra_fix_enforce(LP::Relaxation &rel, int branch_ind,
                        double tour_entry);

/// Undo a Contra branch on an edge which was done by fixing bounds.
void contra_fix_undo(LP::Relaxation &rel, int branch_ind, double tour_entry);

/// Enforce a Contra branch by adding a single inequality.
void contra_naive_enforce(LP::Relaxation &rel, int branch_ind,
                          double tour_entry);

/// Undo a Contra branch on an edge which was done by adding an inequality.
void contra_naive_undo(LP::Relaxation &rel, int branch_ind, double tour_entry);

/// Get the enforcing function that corresponds to \p strat.
std::function<void(LP::Relaxation&, int, double)> enforcer(ContraStrat strat);

/// Get the undoing function that corresponds to \p strat.
std::function<void(LP::Relaxation&, int, double)> undoer(ContraStrat strat);

}
}


#endif
