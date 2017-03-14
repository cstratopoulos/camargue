#include "abc_nodesel.hpp"
#include "graph.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::exception;

using std::vector;

namespace CMR {
namespace ABC {

BoundBrancher::BoundBrancher(const Data::Instance &inst,
                           const LP::ActiveTour &active_tour,
                           const Data::BestGroup &best_data,
                           const Graph::CoreGraph &core_graph,
                           LP::CoreLP &core_lp)
try : BaseBrancher(inst, active_tour, best_data, core_graph, core_lp),
      prob_q(BranchNode::bound_worse)
{
    branch_history.emplace_front();
    prob_q.push(branch_history.begin());
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("BoundBrancher constructor failed.");
}

/**
 * Make a heap insertion of the nodes in \p prob_array to ensure with respect
 * to an ordering that prioritizes unvisited nodes with good tours.
 */
void BoundBrancher::enqueue_split(BranchNode::Split prob_array) try
{
    for (BranchNode &B : prob_array) {
        branch_history.emplace_front(std::move(B));
        prob_q.push(branch_history.begin());
    }
} catch (const exception &e) {
    cerr << e.what() << " putting nodes in history" << endl;
    throw runtime_error("BoundBrancher::enqueue_split failed.");
}

BranchHistory::iterator BoundBrancher::next_prob()
{
    BranchHistory::iterator result = peek_next();
    prob_q.pop();

    return result;
}

BranchHistory::iterator BoundBrancher::peek_next()
{
    if (branch_history.empty() || prob_q.empty())
        return branch_history.end();

    return prob_q.top();
}

}
}
