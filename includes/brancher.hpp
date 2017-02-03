#ifndef CMR_BRANCHER_H
#define CMR_BRANCHER_H

#include "core_lp.hpp"
#include "graph.hpp"
#include "branch_util.hpp"

#include <iostream>
#include <functional>
#include <stack>
#include <vector>

namespace CMR {
namespace ABC {

/// Generation and dispensing of branching subproblems.
class Brancher {
public:
    Brancher(LP::Relaxation &lp_rel, const std::vector<Graph::Edge> &edges,
             const LP::TourBasis &tbase, const double &tourlen,
             const ContraStrat strat);

    /// Pick an edge/branching problem to examine, modifying the Relaxation.
    Problem next_prob();

    /// Mark the top problem as \p stat, undoing and splitting if neccessary.
    void pop_problem(const Problem::Status stat);

    /// Just get the next edge to branch on, without modifying the Relaxation.
    int branch_edge_index();

    bool solved(const Problem prob) const
        { return prob == solved_prob; }

private:

    void split_prob(int edge); //<! Add two child subproblems to the queue.
    void enact_top(); //<! Adjust the lp_relax based on the top node Status.
    
    LP::Relaxation &lp_relax;
    const std::vector<Graph::Edge> &core_edges;
    const LP::TourBasis &tour_base;
    const double &tour_len;

    const ContraStrat contra_strategy;
    const std::function<void(LP::Relaxation&, int, double)> contra_enforce;
    const std::function<void(LP::Relaxation&, int, double)> contra_undo;
    
    std::stack<Problem> subprobs; //<! The problems to be considered.
    static Problem solved_prob;
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
