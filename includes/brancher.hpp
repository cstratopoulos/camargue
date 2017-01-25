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

/// A simple structure for recording the status of branching subproblems.
struct Problem {
    enum class Type {
        Root, //<! The root lp.
        Affirm, //<! Enforcing agreement with current tour.
        Contra //<! Enforcing departure from current tour.
    };

    enum Status {
        Unseen, //<! Unexamined subproblem.
        Seen, //<! Examinded subproblem.
        Pruned //<! Examined and pruned. 
    };
    
    Problem() = default;
    
    Problem(Type ptype, int index)
        : type(ptype), status(Status::Unseen), edge_ind(index) {}
    
    Problem(Type ptype, Status pstat, int index)
        : type(ptype), status(pstat), edge_ind(index) {}

    bool operator==(const Problem &rhs) const
        {
            return (type == rhs.type &&
                    status == rhs.status &&
                    edge_ind == rhs.edge_ind);
        }

    Type type;
    Status status;
    int edge_ind;
};

std::ostream &operator<<(std::ostream &os, Problem::Type type);
std::ostream &operator<<(std::ostream &os, Problem::Status stat);
std::ostream &operator<<(std::ostream &os, const Problem &prob);

/// Strategies for enforcing Problem::Type::Contra branches.
enum class ContraStrat {
    Fix, /// Change the bounds on an edge in the Relaxation.
    Dive, /// Perturb objective function coefficients. 
};

/// Enforce a Contra branch on an edge by fixing bounds.
void contra_fix_enforce(LP::Relaxation &rel, const int branch_ind,
                        const double tour_entry);

/// Undo a Contra branch on an edge which was done by fixing bounds.
void contra_fix_undo(LP::Relaxation &rel, const int branch_ind,
                     const double tour_entry);


/// Generation and dispensing of branching subproblems.
class Brancher {
public:
    Brancher(LP::Relaxation &lp_rel, const std::vector<Graph::Edge> &edges,
             const LP::TourBasis &tbase, const double &tourlen,
             const ContraStrat strat);

    /// Pick an edge/branching problem to examine, modifying the Relaxation.
    Problem &next_prob();

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
    const std::function<void(LP::Relaxation&, const int,
                             const double)> contra_enforce;
    const std::function<void(LP::Relaxation&, const int,
                             const double)> contra_undo;
    
    std::stack<Problem> subprobs; //<! The problems to be considered.
    static Problem solved_prob;
};

}
}


#endif
