#include "exec_branch.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

extern "C" {
#include <concorde/INCLUDE/linkern.h>
#include <concorde/INCLUDE/edgegen.h>
}

using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {
namespace ABC {

Executor::Executor(const Data::Instance &inst,
                   const Data::BestGroup &bestdata,
                   const Graph::CoreGraph &coregraph,
                   LP::CoreLP &core, BranchTourFind &btourfind)
    : instance(inst), active_tour(core.get_active_tour()), best_data(bestdata),
      core_graph(coregraph), core_lp(core), btour_find(btourfind)
{}

/// @returns a ScoreTuple with info about the next edge to add to branch on.
ScoreTuple Executor::branch_edge()
{
    runtime_error err("Problem in Executor::branch_edge");

    vector<double> x = core_lp.lp_vec();
    vector<int> colstat = core_lp.col_stat();
    vector<int> lw_inds;

    for (int i = 0; i < colstat.size(); ++i)
        if (colstat[i] == 1 && !util::var_integral(x[i]))
            lw_inds.push_back(i);

    if (lw_inds.empty())
        throw runtime_error("Tried to branch with no fractional basic vars");

    const vector<Graph::Edge> &core_edges = core_graph.get_edges();
    vector<int> sb1inds;

    try {
        sb1inds = length_weighted_cands(core_edges, lw_inds, x, SB::Cands1);
    } CMR_CATCH_PRINT_THROW("getting longedge candidates", err);

    vector<LP::Estimate> down_ests;
    vector<LP::Estimate> up_ests;
    vector<LP::Basis> cbases;
    double upper_bound = best_data.min_tour_value;
    double avg_itcount = core_lp.avg_itcount();

    try {
        core_lp.primal_strong_branch(active_tour.edges(),
                                     active_tour.base().colstat,
                                     active_tour.base().rowstat,
                                     sb1inds, down_ests, up_ests, cbases,
                                     SB::round1_limit(avg_itcount),
                                     upper_bound);
    } CMR_CATCH_PRINT_THROW("getting first round candidates", err);

    vector<ScoreTuple> sb2cands;
    vector<int> sb2inds;
    vector<LP::Basis> sb2bases;

    try {
        sb2cands = ranked_cands(sb1inds, down_ests, up_ests, core_edges, cbases,
                                SB::StrongMult, upper_bound, SB::Cands2);
    } CMR_CATCH_PRINT_THROW("ranking first round candidates", err);

    try {
        for (ScoreTuple &t : sb2cands) {
            int edge_ind = core_graph.find_edge_ind(t.ends.end[0],
                                                    t.ends.end[1]);
            if (edge_ind == -1)
                throw runtime_error("Candidate edge not in graph");
            sb2inds.push_back(edge_ind);
            sb2bases.emplace_back(std::move(t.contra_base));
        }

        core_lp.primal_strong_branch(active_tour.edges(),
                                     active_tour.base().colstat,
                                     active_tour.base().rowstat,
                                     sb2inds, down_ests, up_ests, sb2bases,
                                     SB::round2_limit(avg_itcount),
                                     upper_bound);
    } CMR_CATCH_PRINT_THROW("getting round 2 cands", err);

    ScoreTuple winner;

    try {
        winner = std::move(ranked_cands(sb2inds, down_ests, up_ests, core_edges,
                                        sb2bases, SB::StrongMult,
                                        upper_bound, 1)[0]);
    } CMR_CATCH_PRINT_THROW("getting winner", err);

    EndPts &wends = winner.ends;
    int winner_ind = core_graph.find_edge_ind(wends.end[0], wends.end[1]);
    if (winner_ind == -1)
        throw runtime_error("Winning branch edge not in graph");

    cout << winner << endl;

    return winner;
}

BranchNode::Split Executor::split_problem(ScoreTuple &branch_tuple,
                                          BranchNode &parent)
{
    using EstStat = LP::Estimate::Stat;

    runtime_error err("Problem in Executor::split_problem");
    const EndPts &branch_edge = branch_tuple.ends;
    vector<EndsDir> edge_stats;

    try {
        edge_stats = btour_find.common_constraints(parent, branch_edge);
    } CMR_CATCH_PRINT_THROW("building edge stats", err);


    BranchNode::Split result;

    for (int i : {0, 1}) {
        edge_stats.emplace_back(EndsDir(branch_edge, dir_from_int(i)));

        bool feas = true;
        double tour_val = 0.0;

        LP::Estimate &est = (i == 0 ? branch_tuple.down_est :
                             branch_tuple.up_est);
        EstStat estat = est.sol_stat;
        double estval = est.value;

        try { btour_find.estimate_tour(edge_stats, feas, tour_val); }
        CMR_CATCH_PRINT_THROW("computing a tour estimate", err);

        result[i] = BranchNode(branch_edge, dir_from_int(i), parent,
                               tour_val, estval);

        if (!feas)
            result[i].stat = BranchNode::Status::Pruned;
        else {
            if (estat != EstStat::Abort || estval > best_data.min_tour_value) {
                result[i].price_basis = std::move(est.sb_base);
                if (estat == EstStat::Infeas)
                    result[i].stat = BranchNode::Status::NeedsRecover;
                else
                    result[i].stat = BranchNode::Status::NeedsPrice;
            }
        }
        edge_stats.pop_back();
    }

    parent.stat = BranchNode::Status::Done;

    return result;
}

void Executor::clamp(const BranchNode &current_node)
{
    if (current_node.is_root())
        return;

    int index = core_graph.find_edge_ind(current_node.ends.end[0],
                                         current_node.ends.end[1]);
    if (index == -1)
        throw runtime_error("Tried to clamp edge not in graph.");

    double clamp_bd = static_cast<int>(current_node.direction);
    char clamp_sense = '\0';

    if (clamp_bd == 0.0)
        clamp_sense = 'U';
    else
        clamp_sense = 'L';

    try {
        core_lp.tighten_bound(index, clamp_sense, clamp_bd);
    } catch (const exception &e) {
        cerr << e.what() << endl;
        throw runtime_error("Executor::clamp failed.");
    }

    if (verbose)
        cout << "Called clamp on " << bnode_brief(current_node)
             << ", set bd " << clamp_sense << " to "
             << static_cast<int>(clamp_bd) << endl;
}

void Executor::unclamp(const BranchNode &current_node)
{
    if (current_node.is_root())
        return;

    int index = core_graph.find_edge_ind(current_node.ends.end[0],
                                         current_node.ends.end[1]);
    if (index == -1)
        throw runtime_error("Tried to unclamp edge not in graph.");

    double clamp_bd = static_cast<int>(current_node.direction);
    char clamp_sense = '\0';

    if (clamp_bd == 0.0)
        clamp_sense = 'U';
    else
        clamp_sense = 'L';

    try {
        core_lp.tighten_bound(index, clamp_sense, 1 - clamp_bd);
    } catch (const exception &e) {
        cerr << e.what() << endl;
        throw runtime_error("Executor::clamp failed.");
    }

    if (verbose)
        cout << "Called unclamp on " << bnode_brief(current_node)
             << ", set bd " << clamp_sense << " to "
             << static_cast<int>(1 - clamp_bd) << endl;
}



}
}
