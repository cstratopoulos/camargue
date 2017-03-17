#include "base_brancher.hpp"

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

BaseBrancher::BaseBrancher(const Data::Instance &inst,
                         const LP::ActiveTour &activetour,
                         const Data::BestGroup &bestdata,
                         const Graph::CoreGraph &coregraph, LP::CoreLP &core)
try : instance(inst), best_data(bestdata), core_graph(coregraph), core_lp(core),
      exec(inst, activetour, bestdata, coregraph, core),
      next_itr(branch_history.end())
{
    branch_history.emplace_front();
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("BaseBrancher constructor failed.");
}

void BaseBrancher::split_prob(BranchHistory::iterator &current)
{
    runtime_error err("Problem in BaseBrancher::next_level");

    BranchNode::Split prob_array;
    ScoreTuple branch_tuple;
    bool verbose = true;

    try {
        cout << "Getting next branch edge...." << endl;
        branch_tuple = exec.branch_edge();
        cout << "Splitting on edge...." << endl;
        prob_array = exec.split_problem(branch_tuple, *current);
    } CMR_CATCH_PRINT_THROW("finding next edge and splitting", err);

    try {
        enqueue_split(std::move(prob_array));
    } CMR_CATCH_PRINT_THROW("adding child problems to queue", err);

    if (verbose) {
        cout << "BaseBrancher::split_prob printing history thus far" << endl;
        for (const auto &B : branch_history)
            cout << B << "\n";
    }
}

void BaseBrancher::do_branch(const BranchNode &B)
{
    if (B.is_root())
        return;

    cout << "Calling do_branch on " << bnode_brief(B) << "..." << endl;

    runtime_error err("Problem in BaseBrancher::do_branch");
    vector<int> tour;

    try {
        const BranchNode *iter = &B;
        vector<EndsDir> edge_stats;

        while(!iter->is_root()) {
            edge_stats.emplace_back(iter->ends, iter->direction);
            iter = iter->parent;
        }

        double tval = 0.0;
        bool found_tour = false;
        bool feas = true;

        exec.branch_tour(edge_stats, core_lp.get_active_tour().nodes(),
                         found_tour, feas,
                         tour, tval, true);
        if (!found_tour)
            throw runtime_error("Unimplemented case of no tour");

    } CMR_CATCH_PRINT_THROW("building edge_stats/getting branch tour", err);

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

void BaseBrancher::do_unbranch(const BranchNode &B)
{
    if (B.is_root())
        return;

    if (!B.visited()) {
        std::string es = "Calling do_unbranch on unvisited " + bnode_brief(B);
        throw runtime_error(es);
    }

    runtime_error err("Problem in BaseBrancher::do_unbranch");

    fetch_next();
    const BranchHistory::iterator Bnext = next_itr;
    if (Bnext == branch_history.end())
        return;

    try {
        common_prep_next(B, *Bnext);
    } CMR_CATCH_PRINT_THROW("calling common_prep_next", err);
}

/**
 * @param done the BranchNode that was just examined/pruned.
 * @param next the BranchNode that will succeed \p done in the current search.
 * Compute the common ancestor A of \p done and \p next, undoing all the clamps
 * from \p done to A and doing all the clamps from A to the parent of
 * \p next.
 */
void BaseBrancher::common_prep_next(const BranchNode &done,
                                   const BranchNode &next)
{
    runtime_error err("Problem in BaseBrancher::common_prep_next");
    cout << "Calling common_prep next.\n The Done node:\n"
         << done << "\nThe Next node:\n" << next << endl;

    if (next.parent == &done) {
        cout << "next is child of done, returning" << endl;
        return;
    }

    try {
        exec.unclamp(done);
    } CMR_CATCH_PRINT_THROW("undoing done problem", err);
    cout << "Undid clamp on done." << endl;

    if (done.parent->is_root()) {
        cout << "done node is root child, common ancestor is root." << endl;
    }

    if (done.parent == next.parent) {
        cout << "done and next are siblings, returning" << endl;
        return;
    }

    using IterPath = std::pair<const BranchNode *, vector<const BranchNode *>>;

    IterPath done_path(done.parent, {});
    IterPath next_path(next.parent, {});

    if (done.parent->depth != next.parent->depth) {
        int shallow_common = std::min(done.parent->depth, next.parent->depth);
        IterPath &catchup = (done.parent->depth > shallow_common ?
                                   done_path : next_path);

        try {
            while (catchup.first->depth > shallow_common &&
                   !catchup.first->is_root()) {
                catchup.second.push_back(catchup.first);
                catchup.first = catchup.first->parent;
            }
        } CMR_CATCH_PRINT_THROW("building catchup list", err);
    }

    try {
        while (done_path.first != next_path.first) {
            done_path.second.push_back(done_path.first);
            next_path.second.push_back(next_path.first);

            done_path.first = done_path.first->parent;
            next_path.first = next_path.first->parent;
        }
    } CMR_CATCH_PRINT_THROW("finding common ancestor path", err);

    cout << "Loop completed with common ancestor "
         << bnode_brief(*(done_path.first)) << endl;

    try {
        cout << "Undoing clamps from done to ancestor..." << endl;
        for (const BranchNode *b : done_path.second) {
            exec.unclamp(*b);
        }

        cout << "Clamping from next to ancestor...." << endl;
        for (const BranchNode *b : next_path.second) {
            exec.clamp(*b);
        }
    } CMR_CATCH_PRINT_THROW("doing actual clamps/unclamps", err);
}

}
}
