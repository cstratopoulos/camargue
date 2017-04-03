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
                             const Data::BestGroup &best_data,
                             const Graph::CoreGraph &core_graph,
                             LP::CoreLP &core_lp)
try : BaseBrancher(inst, best_data, core_graph, core_lp)
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
    cout << "Calling InterBrancher::next_prob...." << endl;
    if (next_itr == branch_history.end()) {
        fetch_next();
    } else {
        cout << "...Already set." << endl;
    }

    BranchHistory::iterator result = next_itr;

    next_itr = branch_history.end();
    return result;
}

void InterBrancher::fetch_next()
{
    using QueueIterator = std::vector<BranchHistory::iterator>::iterator;

    if (verbose)
        cout << "Calling InterBrancher::fetch_next, "
             << prob_q.size() << " problems in q" << endl;
    if (branch_history.empty() || prob_q.empty()) {
        next_itr = branch_history.end();
        return;
    }

    if (node_num % BestFreq != 0) { //pop a usual best tour node from heap
        heap_pop(prob_q);
        next_itr = prob_q.back();
        if (verbose)
            cout << "USUAL (prob " << node_num << ") set next_itr to "
                 << bnode_brief(*next_itr) << endl;
        prob_q.pop_back();
    } else { // pick a best bound node, corrupting and rebuilding the heap.
        QueueIterator best_bound_itr = std::min_element(prob_q.begin(),
                                                        prob_q.end(),
                                                        better_bound);
        QueueIterator last_itr = std::prev(prob_q.end());

        next_itr = *best_bound_itr;
        if (verbose)
            cout << "INTERLEAVE (prob " << node_num << ") set next_itr to "
                 << bnode_brief(*next_itr) << endl;
        std::swap(*best_bound_itr, *last_itr);
        prob_q.pop_back();
        heap_make(prob_q);
    }

    ++node_num;
}

}
}
