#include "dfs_brancher.hpp"
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
try : instance(inst), best_data(bestdata), core_graph(coregraph), core_lp(core),
      exec(inst, activetour, bestdata, coregraph, core)
{
    branch_history.emplace_front();
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("DFSbrancher constructor failed.");
}

SplitIter DFSbrancher::next_level()
{
    runtime_error err("Problem in DFSbrancher::next_level");

    SplitIter result{{branch_history.end(), branch_history.end()}};

    BranchHistory::iterator parent_it = next_prob();
    if(parent_it == branch_history.end()) {
        cout << "No more unvisited nodes, search complete." << endl;
        return result;
    } else {
        cout << "\ngetting parent_it in next level...\n"
             << (*parent_it) << "\n\tvisited: "
             << parent_it->visited() << endl;
    }

    BranchNode::Split prob_array;
    ScoreTuple branch_tuple;

    try {
        cout << "Getting next branch edge...." << endl;
        branch_tuple = exec.branch_edge();
        cout << "Splitting on edge...." << endl;
        prob_array = exec.split_problem(branch_tuple, *parent_it);
    } CMR_CATCH_PRINT_THROW("finding next edge and splitting", err);

    const EndPts branch_edge = branch_tuple.ends;
    int ind = core_graph.find_edge_ind(branch_edge.end[0], branch_edge.end[1]);
    if (ind == -1)
        throw runtime_error("Couldn't get index of branch edge");

    int tour_entry = best_data.best_tour_edges[ind];

    if (tour_entry == 0)
        std::swap(prob_array[0], prob_array[1]);

    try {
        for (BranchNode &B : prob_array)
            branch_history.emplace_front(std::move(B));
    } CMR_CATCH_PRINT_THROW("putting nodes in history", err);

    result[0] = branch_history.begin();
    result[1] = std::next(branch_history.begin());

    cout << "Branch stack thus far:\n";
    for (const auto &B : branch_history)
        cout << B << "\n";


    return result;
}

void DFSbrancher::do_branch(BranchNode &B)
{
    runtime_error err("Problem in DFSbrancher::do_branch");
    vector<int> tour;

    try {
        if (!B.tour_clq)
            throw runtime_error("Null tour clique");

        tour = exec.expand_tour(B.tour_clq);
    } CMR_CATCH_PRINT_THROW("expanding branch tour", err);

    int ncount = core_graph.node_count();
    vector<Graph::Edge> missing_edges;

    try {
        for (int i = 0; i < ncount; ++i) {
            int e0 = tour[i];
            int e1 = tour[(i + 1) % ncount];
            int ind = core_graph.find_edge_ind(e0, e1);
            if (ind == -1)
                missing_edges.emplace_back(e0, e1, instance.edgelen(e0, e1));
        }

        core_lp.add_edges(missing_edges, false);
    } CMR_CATCH_PRINT_THROW("finding add adding missing edges", err);


    try {
        exec.clamp(B);
        core_lp.set_active_tour(std::move(tour));
    } CMR_CATCH_PRINT_THROW("clamping bound/instating branch tour", err);

    cout << "\nBranched " << B << endl;
}

void DFSbrancher::do_unbranch(BranchNode &B)
{
    runtime_error err("Problem in DFSbrancher::do_unbranch");
    vector<int> tour;

    try {
        if (B.parent->is_root())
            tour = best_data.best_tour_nodes;
        else {
            const Sep::Clique::Ptr &parent_tcliq = B.parent->tour_clq;
            if (!parent_tcliq)
                throw runtime_error("Null parent tour");
            tour = exec.expand_tour(parent_tcliq);
        }
    } CMR_CATCH_PRINT_THROW("getting parent tour", err);

    if (core_lp.dual_feas())
        B.stat = BranchNode::Status::Pruned;
    else
        B.stat = BranchNode::Status::Done;

    try {
        exec.unclamp(B);
        core_lp.set_active_tour(std::move(tour));
    } CMR_CATCH_PRINT_THROW("undoing bound/instating parent tour", err);

    cout << "\nUnbranched "
         << bnode_brief(B) << ", "
         << "back to parent tour of length "
         << core_lp.active_tourlen() << endl;
}

BranchHistory::iterator DFSbrancher::next_prob()
{
    return std::find_if_not(branch_history.begin(), branch_history.end(),
                            [](const BranchNode &B)
                            { return B.visited(); });
}

}
}
