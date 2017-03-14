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
                         const LP::ActiveTour &activetour,
                         const Data::BestGroup &bestdata,
                         const Graph::CoreGraph &coregraph, LP::CoreLP &core)
try : BaseBrancher(inst, activetour, bestdata, coregraph, core)
{} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("DFSbrancher constructor failed.");
}

void DFSbrancher::enqueue_split(BranchNode::Split prob_array)
{
    runtime_error err("Problem in DFSbrancher::enqueue_split");

    // The tour entry checking below ensures that the first node examined is
    // always the one that agrees with the current tour.

    const EndPts &branch_edge = prob_array[0].ends;
    int ind = core_graph.find_edge_ind(branch_edge.end[0], branch_edge.end[1]);

    if (ind == -1) {
        cerr << "Branch edge " << branch_edge << " not found in graph" << endl;
        throw err;
    }

    int tour_entry = best_data.best_tour_edges[ind];

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
