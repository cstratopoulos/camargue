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

void InterBrancher::heap_push(vector<BranchHistory::iterator> &target_q,
                              BranchHistory::iterator itr)
{
    target_q.push_back(itr);
    std::push_heap(target_q.begin(), target_q.end(), BranchNode::tour_worse);
}

void InterBrancher::heap_make(vector<BranchHistory::iterator> &target_q)
{
    std::make_heap(target_q.begin(), target_q.end(), BranchNode::tour_worse);
}

void InterBrancher::heap_pop(vector<BranchHistory::iterator> &target_q)
{
    std::pop_heap(target_q.begin(), target_q.end(), BranchNode::tour_worse);
}

InterBrancher::InterBrancher(const Data::Instance &inst,
                             const LP::ActiveTour &active_tour,
                             const Data::BestGroup &best_data,
                             const Graph::CoreGraph &core_graph,
                             LP::CoreLP &core_lp)
try : BaseBrancher(inst, active_tour, best_data, core_graph, core_lp)
{
    heap_push(prob_q, branch_history.begin());
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("InterBrancher constructor failed.");
}

/**
 * Make a heap insertion of the nodes in \p prob_array to ensure with respect
 * to an ordering that prioritizes unvisited nodes with good tours.
 */
void InterBrancher::enqueue_split(BranchNode::Split prob_array) try
{
    for (BranchNode &B : prob_array) {
        branch_history.emplace_front(std::move(B));
        heap_push(prob_q, branch_history.begin());
    }
} catch (const exception &e) {
    cerr << e.what() << " putting nodes in history" << endl;
    throw runtime_error("InterBrancher::enqueue_split failed.");
}

BranchHistory::iterator InterBrancher::next_prob()
{
    BranchHistory::iterator result = peek_next();
    prob_q.pop_back();

    if (node_num % BestFreq == 0) { //then the heap needs to be rebuilt.
        if (result != branch_history.end())
            cout << "INTERLEAVE: returning " << bnode_brief(*result) << endl;
        heap_make(prob_q);
    } else if (result != branch_history.end())
        cout << "USUAL (prob " << node_num << "):"
             << " returning " << bnode_brief(*result) << endl;

    ++node_num;
    return result;
}

BranchHistory::iterator InterBrancher::peek_next()
{
    if (branch_history.empty() || prob_q.empty())
        return branch_history.end();

    if (node_num % BestFreq != 0) { //pick a usual best tour node.
        heap_pop(prob_q);
    } else { // pick a best bound node, corrupting the heap.
        auto best_bound_itr = std::min_element(prob_q.begin(),
                                               prob_q.end(), better_bound);
        auto last_itr = std::prev(prob_q.end());

        if (best_bound_itr != last_itr)
            std::swap(*best_bound_itr, *last_itr);
    }

    return prob_q.back();
}

}
}
