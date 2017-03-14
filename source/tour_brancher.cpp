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

TourBrancher::TourBrancher(const Data::Instance &inst,
                           const LP::ActiveTour &active_tour,
                           const Data::BestGroup &best_data,
                           const Graph::CoreGraph &core_graph,
                           LP::CoreLP &core_lp)
try : BaseBrancher(inst, active_tour, best_data, core_graph, core_lp)
      {} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("TourBrancher constructor failed.");
}

/**
 * Make a heap insertion of the nodes in \p prob_array to ensure with respect
 * to an ordering that prioritizes unvisited nodes with good tours.
 */
void TourBrancher::enqueue_split(BranchNode::Split prob_array) try
{
    for (BranchNode &B : prob_array) {
        branch_history.emplace_front(std::move(B));
        // std::push_heap(branch_history.begin(), branch_history.end(),
        //                BranchNode::tour_compare);
    }
} catch (const exception &e) {
    cerr << e.what() << " putting nodes in history" << endl;
    throw runtime_error("TourBrancher::enqueue_split failed.");
}

BranchHistory::iterator TourBrancher::next_prob()
{
    if (branch_history.empty())
        return branch_history.end();

    cout << "Calling TourBrancher::next_prob...." << endl;

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

    if (result->visited()) {
        cout << "All problems visited, returning end." << endl;
        return branch_history.end();
    } else {
        cout << "Returning " << bnode_brief(*result) << ", tour length "
             << static_cast<int>(result->tourlen) << endl;
        return result;
    }
}

}
}
