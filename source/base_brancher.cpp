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
                           const Data::BestGroup &bestdata,
                           const Graph::CoreGraph &coregraph, LP::CoreLP &core)
try : instance(inst), best_data(bestdata), core_graph(coregraph), core_lp(core),
      btour_find(inst, bestdata, coregraph, core),
      exec(inst, bestdata, coregraph, core, btour_find),
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

    try {
        if (verbose)
            cout << "Getting next branch edge...." << endl;
        branch_tuple = exec.branch_edge();
        if (verbose)
            cout << "Splitting on edge...." << endl;
        prob_array = exec.split_problem(branch_tuple, *current);
    } CMR_CATCH_PRINT_THROW("finding next edge and splitting", err);

    try {
        enqueue_split(std::move(prob_array));
    } CMR_CATCH_PRINT_THROW("adding child problems to queue", err);

    int num_remain = 0;

    if (verbose)
        cout << "BaseBrancher::split_prob printing history thus far"
             << endl;

    for (const auto &B : branch_history) {
        if (verbose)
            cout << B << "\n";
        if (!B.visited())
            ++num_remain;
    }

    cout << num_remain << " unvisited problems remain" << endl;
}

void BaseBrancher::do_branch(const BranchNode &B)
{
    if (B.is_root())
        return;

    if (verbose)
        cout << "Calling do_branch on " << bnode_brief(B) << "..." << endl;

    runtime_error err("Problem in BaseBrancher::do_branch");
    vector<int> tour;
    bool found_tour = false;

    try { btour_find.instate_branch_tour(B, found_tour, tour); }
    CMR_CATCH_PRINT_THROW("getting branch tour", err);

    try {
        exec.clamp(B);
        if (found_tour)
            core_lp.set_active_tour(std::move(tour));
    } CMR_CATCH_PRINT_THROW("clamping bound/instating branch tour", err);

    if (verbose)
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

    if (verbose > 1)
        cout << "Calling common_prep next.\n The Done node:\n"
             << done << "\nThe Next node:\n" << next << endl;

    if (next.parent == &done) {
        if (verbose > 1)
            cout << "next is child of done, returning" << endl;
        return;
    }

    try {
        exec.unclamp(done);
    } CMR_CATCH_PRINT_THROW("undoing done problem", err);

    if (verbose > 1)
        cout << "Undid clamp on done." << endl;

    if (done.parent->is_root()) {
        if (verbose > 1)
            cout << "done node is root child, common ancestor is root." << endl;
    }

    if (done.parent == next.parent) {
        if (verbose > 1)
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

    if (verbose > 1)
        cout << "Loop completed with common ancestor "
             << bnode_brief(*(done_path.first)) << endl;

    try {
        if (verbose > 1)
            cout << "Undoing clamps from done to ancestor..." << endl;
        for (const BranchNode *b : done_path.second) {
            exec.unclamp(*b);
        }

        if (verbose > 1)
            cout << "Clamping from next to ancestor...." << endl;
        for (const BranchNode *b : next_path.second) {
            exec.clamp(*b);
        }
    } CMR_CATCH_PRINT_THROW("doing actual clamps/unclamps", err);
}

}
}
