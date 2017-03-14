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

DFSbrancher::DFSbrancher(const Data::Instance &inst,
                         const LP::ActiveTour &active_tour,
                         const Data::BestGroup &best_data,
                         const Graph::CoreGraph &core_graph,
                         LP::CoreLP &core_lp)
try : BaseBrancher(inst, active_tour, best_data, core_graph, core_lp)
{} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("DFSbrancher constructor failed.");
}

/**
 * If applicable, this method will swap the position of nodes in
 * \p prob_array to ensure that the node agreeing with the current tour is
 * first in the queue.
 */
void DFSbrancher::enqueue_split(BranchNode::Split prob_array)
{
    runtime_error err("Problem in DFSbrancher::enqueue_split");

    const EndPts &branch_edge = prob_array[0].ends;
    int ind = core_graph.find_edge_ind(branch_edge.end[0], branch_edge.end[1]);

    if (ind == -1) {
        cerr << "Branch edge " << branch_edge << " not found in graph" << endl;
        throw err;
    }

    int tour_entry = best_data.best_tour_edges[ind];

    // swap based on tour entry since they will be added to the list in
    // the reverse order of their appearance in the array.
    if (tour_entry == 0) {
        std::swap(prob_array[0], prob_array[1]);
    }

    try {
        for (BranchNode &B : prob_array)
            branch_history.emplace_front(std::move(B));
    } CMR_CATCH_PRINT_THROW("putting nodes in history", err);
}


BranchHistory::iterator DFSbrancher::next_prob()
{
    return std::find_if_not(branch_history.begin(), branch_history.end(),
                            [](const BranchNode &B)
                            { return B.visited(); });
}

}
}
