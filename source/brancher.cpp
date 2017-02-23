#include "brancher.hpp"
#include "err_util.hpp"
#include "util.hpp"

extern "C" {
#include <concorde/INCLUDE/linkern.h>
}

#include <memory>
#include <stdexcept>

using std::function;
using std::unique_ptr;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {
namespace ABC {

using Edge = Graph::Edge;

using Ptype = Problem::Type;

using Strat = ContraStrat;

Brancher::Brancher(LP::Relaxation &lp_rel,
                   const Graph::CoreGraph &coregraph,
                   const Data::BestGroup &bestdata,
                   const LP::ActiveTour &activetour,
                   const ContraStrat strat) try
    : lp_relax(lp_rel),
      core_graph(coregraph), best_data(bestdata),
      active_tour(activetour),
      contra_strategy(strat),
      contra_enforce(enforcer(strat)), contra_undo(undoer(strat))
{} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Brancher constructor failed.");
}

std::array<Problem, 2> Brancher::next_level()
{
    runtime_error err("Problem in Brancher::next_level");

    ScoreTuple next_obj;

    try { next_obj = next_branch_obj(); }
    CMR_CATCH_PRINT_THROW("getting branching structure", err);

    int index = next_obj.index;
    double tour_entry = best_data.best_tour_edges[index];

    if (tour_entry == 0) {
        return std::array<Problem, 2>{Problem(index, next_obj.down_est),
            Problem(index, next_obj.up_est, std::move(next_obj.contra_base))};
    } else {
        return std::array<Problem, 2>{Problem(index, next_obj.up_est),
            Problem(index, next_obj.down_est, std::move(next_obj.contra_base))};
    }
}

ScoreTuple Brancher::next_branch_obj()
{
    vector<double> x = lp_relax.lp_vec();
    vector<int> colstat = lp_relax.col_stat();
    vector<int> lw_inds;

    for (int i = 0; i < colstat.size(); ++i)
        if (colstat[i] == 1 && !util::var_integral(x[i]))
            lw_inds.push_back(i);

    if (lw_inds.empty())
        throw logic_error("Tried to branch with no fractional basic vars");

    const vector<Graph::Edge> &core_edges = core_graph.get_edges();
    vector<int> sb1inds = length_weighted_cands(core_edges, lw_inds, x,
                                                SB1Cands);
    vector<ScorePair> downobj;
    vector<ScorePair> upobj;
    vector<LP::Basis> cbases;
    double tour_len = best_data.min_tour_value;

    lp_relax.primal_strong_branch(active_tour.edges(),
                                  active_tour.base().colstat,
                                  active_tour.base().rowstat,
                                  sb1inds, downobj, upobj, cbases,
                                  SB1Lim, tour_len);

    vector<ScoreTuple> sb2cands = ranked_cands(sb1inds, downobj, upobj, cbases,
                                               StrongMult, tour_len, SB2Cands);
    vector<int> sb2inds;
    vector<LP::Basis> sb2bases;

    for (ScoreTuple &t : sb2cands) {
        sb2inds.push_back(t.index);
        sb2bases.emplace_back(std::move(t.contra_base));
    }

    lp_relax.primal_strong_branch(active_tour.edges(),
                                  active_tour.base().colstat,
                                  active_tour.base().rowstat,
                                  sb2inds, downobj, upobj, sb2bases,
                                  SB2Lim, tour_len);

    ScoreTuple winner = std::move(ranked_cands(sb2inds, downobj, upobj,
                                               sb2bases,
                                               StrongMult, tour_len, 1)[0]);

    cout << "\n\tWinner edge " << winner.index << ", tour entry "
         << best_data.best_tour_edges[winner.index] << "\n";
    cout << "\t\tDown priority " << winner.down_est.first << ", estimate "
         << winner.down_est.second << "\n"
         << "\t\tUp winner priority " << winner.up_est.first << ", estimate "
         << winner.up_est.second << "\n";

    return winner;
}

vector<int> Brancher::branch_tour(const Data::Instance &inst,
                                  vector<int> &start_tour_nodes)
{
    runtime_error err("Problem in Brancher::branch_tour");

    int fixcount = 0;
    vector<int> fixed_edges;
    vector<Graph::Edge> avoid_edges;

    double bt_time = util::zeit();

    try {
        for (const EdgeStats &stat : branch_stats)
            if (stat.second == 0)
                avoid_edges.push_back(stat.first);
            else {
                ++fixcount;
                fixed_edges.push_back(stat.first.end[0]);
                fixed_edges.push_back(stat.first.end[1]);
            }
    } CMR_CATCH_PRINT_THROW("making record of up/down edges", err);

    cout << "\n\t" << fixcount << " fixed to 1, " << avoid_edges.size()
         << " to avoid." << endl;

    int ncount = inst.node_count();
    vector<Graph::Edge> edges_copy;

    try {
        vector<int> elist;
        vector<int> ecap;

        edges_copy = core_graph.get_edges();
        Graph::get_elist(edges_copy, elist, ecap);

        int default_len = large_len(ncount, ecap);

        for (int i = 0; i < fixcount; ++i) {
            int ind = core_graph.find_edge_ind(fixed_edges[2 * i],
                                               fixed_edges[(2 * i) + 1]);
            if (ind == -1)
                throw logic_error("Fix edge not in graph");

            edges_copy[ind].len = -default_len;
        }

        for (Graph::Edge &e : avoid_edges) {
            int ind = core_graph.find_edge_ind(e.end[0], e.end[1]);
            if (ind == -1)
                throw logic_error("Branched edge not in core");
            edges_copy[ind].removable = true;
        }

        edges_copy.erase(std::remove_if(edges_copy.begin(), edges_copy.end(),
                                        [](const Edge &e)
                                        { return e.removable; }),
                         edges_copy.end());
    } CMR_CATCH_PRINT_THROW("copying edges/making temp Instance", err);

    vector<int> elist;
    vector<int> ecap;
    Data::Instance sparse_inst;

    try {
        Graph::get_elist(edges_copy, elist, ecap);
        sparse_inst = Data::Instance(inst.problem_name(), inst.seed(),
                                     inst.node_count(), elist, ecap);
    } CMR_CATCH_PRINT_THROW("getting actual sparse inst", err);


    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);

    double val = 0;
    int kicks = (ncount > 1000 ? 500 : ncount / 2);

    vector<int> result;
    try { result.resize(ncount); }
    CMR_CATCH_PRINT_THROW("allocating result", err)

    cout << "\tCalling linkern tour...." << endl;
    if (CClinkern_tour(ncount, sparse_inst.ptr(), ecap.size(),
                       &elist[0], ncount, kicks,
                       &start_tour_nodes[0], &result[0], &val,
                       1, 0, 0,
                       (char *) NULL, CC_LK_RANDOM_KICK, &rstate))
        throw runtime_error("CClinkern_tour failed.");

    double tourval = 0.0;
    int fc_found = 0;

    for (int i = 0; i < ncount; ++i) {
        EndPts e(result[i], result[(i + 1) % ncount]);
        tourval += inst.edgelen(e.end[0], e.end[1]);
        for (Graph::Edge &ae : avoid_edges) {
            if (ae == e) {
                cerr << "Avoided edge " << e << " is in tour" << endl;
                throw err;
            }
        }

        for (int i = 0; i < fixcount; ++i) {
            EndPts fe(fixed_edges[2 * i], fixed_edges[(2 * i) + 1]);
            if (e == fe) {
                ++fc_found;
                break;
            }
        }
    }

    if (fc_found != fixcount) {
        cerr << "Found " << fc_found << " fixed edges, supposed to be "
             << fixcount << endl;
        throw err;
    }

    bt_time = util::zeit() - bt_time;
    cout << "\tComputed branch tour of length " << tourval << " in "
         << bt_time << "s\n" << endl;

    return result;
}

void Brancher::do_branch(Problem &prob)
{
    runtime_error err("Problem in Brancher::do_branch");

    int ind = prob.edge_ind;
    double tour_entry = best_data.best_tour_edges[ind];
    Ptype ptype = prob.type;

    const vector<Graph::Edge> &core_edges = core_graph.get_edges();

    cout << "\tDoing " << prob << ", ";

    if (ptype == Ptype::Affirm) {
        try {
            if (tour_entry == 1) {
                cout << "clamping to 1";
                lp_relax.tighten_bound(ind, 'L', 1);
            } else if (tour_entry == 0) {
                cout << "clamping to 0";
                lp_relax.tighten_bound(ind, 'U', 0);
            }

            branch_stats.emplace_back(EdgeStats(core_edges[ind], tour_entry));
        } CMR_CATCH_PRINT_THROW("doing affirm branch", err);
    } else if (ptype == Ptype::Contra) {
        cout << "calling contra_enforce";
        try {
            contra_enforce(lp_relax, ind, tour_entry);
            branch_stats.emplace_back(EdgeStats(core_edges[ind],
                                                1 - tour_entry));
        } CMR_CATCH_PRINT_THROW("doing contra branch", err);
    }
    cout << "\n\n";
}

void Brancher::undo_branch(Problem &prob)
{
    runtime_error err("Problem in Brancher::undo_branch");

    cout << "\tUndoing " << prob << ", ";

    int ind = prob.edge_ind;
    double tour_entry = best_data.best_tour_edges[ind];
    Ptype ptype = prob.type;

    if (ptype == Ptype::Affirm) {
        try {
            if (tour_entry == 1) {
                cout << "resetting lb to zero";
                lp_relax.tighten_bound(ind, 'L', 0);
            } else if (tour_entry == 0) {
                cout << "resetting ub to one";
                lp_relax.tighten_bound(ind, 'U', 1);
            }
        } CMR_CATCH_PRINT_THROW("undoing affirm branch", err);
    } else if (ptype == Ptype::Contra) {
        cout << "calling contra_undo";
        try {
            contra_undo(lp_relax, ind, tour_entry);
        } CMR_CATCH_PRINT_THROW("undoing contra branch", err);
    }

    const vector<Graph::Edge> &core_edges = core_graph.get_edges();

    branch_stats.erase(std::remove_if(branch_stats.begin(),
                                         branch_stats.end(),
                                         [ind, &core_edges]
                                         (const EdgeStats &e)
                                         {
                                             return e.first.end ==
                                             core_edges[ind].end;
                                         }),
                          branch_stats.end());


    cout << "\n\n";
}


/**
 * @param[in] rel the Relaxation to modify.
 * @param[in] branch_ind the index of the edge to branch on.
 * @param[in] tour_entry the entry of the best tour corresponding to
 * \p branch_ind.
 * This function modifies the Relaxation to explore a branch where the solution
 * is fixed to disagree with the current best tour. This is done by changing
 * bounds on \p branch_ind so that the upper bound is zero if \p tour_entry
 * is one, and the lower bound is one if \p tour_entry is zero. In either case,
 * effectively adds the constraint that `x[branch_ind] = 1 - tour_entry` to
 * \p rel.
 */
void contra_fix_enforce(LP::Relaxation &rel, const int branch_ind,
                             const double tour_entry)
{
    if (tour_entry == 1)
        rel.tighten_bound(branch_ind, 'U', 0);
    else if (tour_entry == 0)
        rel.tighten_bound(branch_ind, 'L', 1);

    return;
}

/**
 * @param[in] rel the Relaxation to modify.
 * @param[in] branch_ind the index of the edge to undo a branch on.
 * @param[in] tour_entry the entry of the best tour corresponding to
 * \p branch_ind.
 * This function undoes a branch which was performed by contra_fix_enforce. It
 * assumes that upper or lower bounds have been modified precisely in the way
 * indicated by that function, hence restores them to normal values. If
 * \p tour_entry is one, then branch_ind will have its upper bound relaxed
 * back to one. If it is zero, it will have its lower bound relaxed back to
 * zero.
 */
void contra_fix_undo(LP::Relaxation &rel, const int branch_ind,
                             const double tour_entry)
{
    if (tour_entry == 1)
        rel.tighten_bound(branch_ind, 'U', 1);
    else if (tour_entry == 0)
        rel.tighten_bound(branch_ind, 'L', 0);

    return;
}

function<void(LP::Relaxation&, int, double)> enforcer(Strat strat)
{
    switch (strat) {
    case Strat::Fix:
        return contra_fix_enforce;
    default:
        throw logic_error("Tried to get unimplemented enforcer.");
    }
}

function<void(LP::Relaxation&, int, double)> undoer(Strat strat)
{
    switch (strat) {
    case Strat::Fix:
        return contra_fix_undo;
    default:
        throw logic_error("Tried to get unimplemented undoer.");
    }
}

}
}
