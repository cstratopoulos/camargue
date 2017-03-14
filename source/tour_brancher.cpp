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

BESbrancher::BESbrancher(const Data::Instance &inst,
                         const LP::ActiveTour &active_tour,
                         const Data::BestGroup &best_data,
                         const Graph::CoreGraph &core_graph,
                         LP::CoreLP &core_lp)
try : BaseBrancher(inst, active_tour, best_data, core_graph, core_lp)
{} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("BESbrancher constructor failed.");
}

/**
 * Make a heap insertion of the nodes in \p prob_array to ensure with respect
 * to an ordering that prioritizes unvisited nodes with good tours.
 */
void BESbrancher::enqueue_split(BranchNode::Split prob_array) try
{
    for (BranchNode &B : prob_array) {
        branch_history.emplace_back(std::move(B));
        // std::push_heap(branch_history.begin(), branch_history.end(),
        //                BranchNode::tour_compare);
    }
} catch (const exception &e) {
    cerr << e.what() << " putting nodes in history" << endl;
    throw runtime_error("BESbrancher::enqueue_split failed.");
}

BranchHistory::iterator BESbrancher::next_prob()
{
    if (branch_history.empty())
        return branch_history.end();

    // std::make_heap(branch_history.begin(), branch_history.end(),
    //                BranchNode::tour_compare);
    // std::pop_heap(branch_history.begin(), branch_history.end(),
    //               BranchNode::tour_compare);

    BranchHistory::iterator result =
    std::min_element(branch_history.begin(),
                     branch_history.end(),
                     [](const BranchNode &A, const BranchNode &B)
                     {
                         return (std::make_tuple(A.visited(),A.tourlen) <
                                 std::make_tuple(B.visited(), B.tourlen));
                     });

    //BranchHistory::iterator result = std::prev(branch_history.end());

    if (result->visited())
        return branch_history.end();
    else
        return result;
}

}
}
