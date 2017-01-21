#ifndef CMR_BRANCHER_H
#define CMR_BRANCHER_H

#include "core_lp.hpp"
#include "graph.hpp"
#include "branch_util.hpp"

#include <iostream>
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

    Type type;
    Status status;
    int edge_ind;
};

std::ostream &operator<<(std::ostream &os, Problem::Type type);
std::ostream &operator<<(std::ostream &os, Problem::Status stat);
std::ostream &operator<<(std::ostream &os, Problem &prob);



/// Generation and dispensing of branching subproblems.
class Brancher {
public:
    Brancher(LP::Relaxation &lp_rel, const std::vector<Graph::Edge> &edges,
             const LP::TourBasis &tbase, const double &tourlen);

    Problem &next_prob();

private:
    int branch_edge(); //<! Get the next edge to branch on.
    void split_prob(int edge); //<! Add two child subproblems to the queue.
    void enact_top(); //<! Adjust the lp_relax based on the top node Status.
    
    LP::Relaxation &lp_relax;
    const std::vector<Graph::Edge> &core_edges;
    const LP::TourBasis &tour_base;
    const double &tour_len;
    
    std::stack<Problem> subprobs; //<! The problems to be considered.
    static Problem solved_prob;
};

}
}


#endif
