#include "exec_branch.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

extern "C" {
#include <concorde/INCLUDE/linkern.h>
#include <concorde/INCLUDE/edgegen.h>
}

using std::vector;

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

constexpr int IntMax = std::numeric_limits<int>::max();

namespace ABC {

BranchNode::BranchNode() : stat(Status::NeedsBranch), parent(nullptr),
                           depth(0),
                           tour_clq(nullptr), tourlen(IntMax),
                           maybe_infeas(false)
{}

BranchNode::BranchNode(EndPts ends_, Dir direction_,
                       const BranchNode &parent_, Sep::Clique::Ptr tour_clq_,
                       double tourlen_)
    : ends(ends_), direction(direction_), stat(Status::NeedsCut),
      parent(&parent_), depth(1 + parent_.depth),
      tour_clq(tour_clq_), tourlen(tourlen_), maybe_infeas(false) {}

BranchNode::BranchNode(BranchNode &&B) noexcept
    : ends(std::move(B.ends)),
      direction(B.direction),
      stat(B.stat),
      parent(B.parent),
      depth(B.depth),
      tour_clq(std::move(B.tour_clq)),
      tourlen(B.tourlen),
      maybe_infeas(B.maybe_infeas)
{
    B.stat = Status::Done;
    B.parent = nullptr;
    B.depth = 0;
    B.tourlen = IntMax;
    B.maybe_infeas = false;
}

BranchNode &BranchNode::operator=(BranchNode &&B) noexcept
{
    ends = std::move(B.ends);
    direction = B.direction;
    stat = B.stat;
    parent = B.parent;
    depth = B.depth;
    tour_clq = std::move(B.tour_clq);
    tourlen = B.tourlen;
    maybe_infeas = B.maybe_infeas;

    B.stat = Status::Done;
    B.parent = nullptr;
    B.tourlen = IntMax;
    B.maybe_infeas = false;

    return *this;
}

BranchNode::Dir dir_from_int(int entry)
{
    switch (entry) {
    case 0:
        return BranchNode::Dir::Down;
    case 1:
        return BranchNode::Dir::Up;
    default:
        throw runtime_error("Tried to get BranchNode::Dir w non-binary int");
    }
}

std::string bnode_brief(const BranchNode &B)
{
    std::stringstream result;

    if (B.is_root())
        result << "[root]";
    else {
        result << "[" << B.ends << " = " << B.direction << "]";
    }

    return result.str();
}

ostream &operator<<(ostream &os, const BranchNode::Dir &dir)
{
    os << static_cast<int>(dir);
    return os;
}
ostream &operator<<(ostream &os, const BranchNode &B)
{
    os << bnode_brief(B) << ", ";

    if (B.is_root()) {
        os << "depth 0";
    } else {
        os << "depth " << B.depth << ", tourlen " << B.tourlen
           << ", parent "
           << bnode_brief(*(B.parent));
    }
    return os;
}

Executor::Executor(const Data::Instance &inst,
                   const LP::ActiveTour &activetour,
                   const Data::BestGroup &bestdata,
                   const Graph::CoreGraph &coregraph,
                   LP::CoreLP &core) try
    : instance(inst), active_tour(activetour), best_data(bestdata),
      core_graph(coregraph), core_lp(core),
      inds_table(inst.node_count()),
      fix_degrees(inst.node_count()), avail_degrees(inst.node_count()),
      tour_cliques(bestdata.best_tour_nodes, bestdata.perm)
{} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Executor constructor failed");
}

ScoreTuple Executor::branch_edge()
{
    runtime_error err("Problem in Executor::branch_edge");
    using LP::InfeasObj;

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

    vector<InfeasObj> downobj;
    vector<InfeasObj> upobj;
    vector<LP::Basis> cbases;
    double upper_bound = best_data.min_tour_value;
    double avg_itcount = core_lp.avg_itcount();

    try {
        core_lp.primal_strong_branch(active_tour.edges(),
                                     active_tour.base().colstat,
                                     active_tour.base().rowstat,
                                     sb1inds, downobj, upobj, cbases,
                                     SB::round1_limit(avg_itcount),
                                     upper_bound);
    } CMR_CATCH_PRINT_THROW("getting first round candidates", err);

    vector<ScoreTuple> sb2cands;
    vector<int> sb2inds;
    vector<LP::Basis> sb2bases;

    try {
        sb2cands = ranked_cands(sb1inds, downobj, upobj, core_edges, cbases,
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
                                     sb2inds, downobj, upobj, sb2bases,
                                     SB::round2_limit(avg_itcount),
                                     upper_bound);
    } CMR_CATCH_PRINT_THROW("getting round 2 cands", err);

    ScoreTuple winner;

    try {
        winner = std::move(ranked_cands(sb2inds, downobj, upobj, core_edges,
                                        sb2bases, SB::StrongMult,
                                        upper_bound, 1)[0]);
    } CMR_CATCH_PRINT_THROW("getting winner", err);

    EndPts &wends = winner.ends;
    int winner_ind = core_graph.find_edge_ind(wends.end[0], wends.end[1]);
    if (winner_ind == -1)
        throw runtime_error("Winning branch edge not in graph");

    cout << "\n\t" << winner << endl;
    cout << "\t\t Tour entry " << best_data.best_tour_edges[winner_ind]
         << endl;
    if (best_data.best_tour_edges[winner_ind] > 1) {
        cout << "tour edges size "
             << best_data.best_tour_edges.size() << endl;
        cout << "core ecount " << core_graph.edge_count() << endl;
        cout << "winner ind " << winner_ind << endl;
        throw logic_error("winner ind out of range or nb tour entry?????");
    }

    return winner;
}

BranchNode::Split Executor::split_problem(const ScoreTuple &branch_tuple,
                                          BranchNode &parent)
{
    runtime_error err("Problem in Executor::split_problem");
    const EndPts &branch_edge = branch_tuple.ends;
    vector<EndsDir> edge_stats;

    const BranchNode *iter = &parent;
    cout << "Building edge stats starting from parent...." << endl;

    try {
        while(!iter->is_root()) {

            cout << bnode_brief(*iter) << ", ";

            if (iter->ends == branch_edge)
                throw runtime_error("Trying to split already-branched edge.");
            edge_stats.push_back(EndsDir(iter->ends, iter->direction));
            iter = iter->parent;
        }
        edge_stats.reserve(edge_stats.size() + 1);
    } CMR_CATCH_PRINT_THROW("building edge_stats", err);

    cout << endl;

    BranchNode::Split result;
    const vector<int> &start_nodes = active_tour.nodes();

    for (int i : {0, 1}) {
        cout << "Adding " << branch_edge << " = " << i << ", "
             << "getting tour..." << endl;

        edge_stats.emplace_back(EndsDir(branch_edge, dir_from_int(i)));

        bool found_tour = true;
        bool feas = true;
        double tour_val = 0.0;
        vector<int> btour;
        Sep::Clique::Ptr tcliq(nullptr);

        try {
            branch_tour(edge_stats, start_nodes,
                        found_tour, feas,
                        btour, tour_val);
            if (found_tour)
                tcliq = compress_tour(btour);
        } CMR_CATCH_PRINT_THROW("computing/compressing branch tour", err);

        result[i] = BranchNode(branch_edge, dir_from_int(i), parent,
                               tcliq, tour_val);

        if (!feas)
            result[i].stat = BranchNode::Status::Pruned;

        edge_stats.pop_back();
    }

    bool down_infeas{branch_tuple.down_est.first};
    bool up_infeas{branch_tuple.up_est.first};

    if (down_infeas && up_infeas) {
        cerr << "Both subproblems infeasible" << endl;
        throw err;
    }

    result[0].maybe_infeas = down_infeas;
    result[1].maybe_infeas = up_infeas;

    parent.stat = BranchNode::Status::Done;

    return result;
}

/**
 * @param[in] edge_stats a vector of edges being branched upon, together with
 * the value they are clamped to.
 * @param[in] start_tour_nodes the tour being used as the starting cycle for
 * computing \p tour.
 * @param[out] found_tour set to true iff a compliant tour was found.
 * @param[out] feas set to false iff the branch instructions were proven
 * infeasible.
 * @param[out] tour the branching tour.
 * @param[out] tour_val the length of \p tour.
 * @result Uses a short run of chained Lin-Kernighan to compute a tour which
 * complies with the branching directions specified in \p edge_stats, with
 * \p start_tour_nodes as the starting cycle for the LK run. The resulting
 * tour can be used as a starting basis for cutting and pivoting on a branching
 * problem which complies with \p edge_stats. The computed tour will be
 * verified for compliance before it is returned. If no compliant tour is found,
 * \p found_tour will be set to false. If it is proved that \p edge_stats make
 * it impossible for any tour to be found, \p feas will be set to false, as
 * will \p found_tour. If either \p feas or \p found_tour is set to false,
 * \p tour will be cleared and \p tour_val will be set to a large positive
 * value.
 */
void Executor::branch_tour(const vector<EndsDir> &edge_stats,
                           const vector<int> &start_tour_nodes,
                           bool &found_tour, bool &feas,
                           vector<int> &tour, double &tour_val)
{
    runtime_error err("Problem in Executor::branch_tour");

    found_tour = true;
    feas = true;

    auto tour_guard = util::make_guard([&found_tour, &feas,
                                        &tour, &tour_val]() -> void
                                       {
                                           if (!found_tour || !feas) {
                                               found_tour = false;
                                               tour.clear();
                                               tour_val = IntMax;
                                           }
                                       });

    int ncount = core_graph.node_count();

    bool found_contra = false;

    vector<int> want_inds;
    vector<int> avoid_inds;

    std::fill(fix_degrees.begin(), fix_degrees.end(), 0);
    std::fill(avail_degrees.begin(), avail_degrees.end(), ncount - 1);

    inds_table.fill('\0');

    try {
        for (const EndsDir &ed : edge_stats) {
            const EndPts &e = ed.first;
            int ind = core_graph.find_edge_ind(e.end[0], e.end[1]);

            if (ind == -1)
                throw runtime_error("Edge not in graph for branch_tour");

            int target_entry = static_cast<int>(ed.second);
            if (best_data.best_tour_edges[ind] != target_entry) {
                found_contra = true;
            }

            if (target_entry == 0) {
                inds_table(e.end[0], e.end[1]) = 'A';
                avoid_inds.push_back(ind);

                for (int pt : e.end) {
                    --avail_degrees[pt];
                    if (avail_degrees[pt] < 2) {
                        cout << pt
                             << "has degree less than 2 in branch graph, "
                             << "tour infeasible branch." << endl;
                        feas = false;
                        return;
                    }
                }
            } else {
                inds_table(e.end[0], e.end[1]) = 'W';
                want_inds.push_back(ind);

                for (int pt : e.end) {
                    ++fix_degrees[pt];
                    if (fix_degrees[pt] > 2) {
                        cout << fix_degrees[pt] << " fixed vxs incident with "
                             << pt << ", tour infeasible branch." << endl;
                        feas = false;
                        return;
                    }
                }
            }
        }
    } CMR_CATCH_PRINT_THROW("categorizing edges", err);

    if (!found_contra) {
        try {
            cout << "\tFixed edges affirm best tour, returning best tour."
                 << endl;
            tour = best_data.best_tour_nodes;
            tour_val = best_data.min_tour_value;
        } CMR_CATCH_PRINT_THROW("copying best tour", err);

        return;
    }

    cout << "\t" << want_inds.size() << " fixed to 1, "
         << avoid_inds.size() << " to avoid." << endl;

    vector<Graph::Edge> edges_copy;

    try {
        edges_copy = core_graph.get_edges();

        int default_len = large_len(ncount, edges_copy);

        for (int ind : want_inds)
            edges_copy[ind].len = -default_len;

        for (int ind : avoid_inds)
            edges_copy[ind].removable = true;

        edges_copy.erase(std::remove_if(edges_copy.begin(), edges_copy.end(),
                                        [](const Graph::Edge &e)
                                        { return e.removable; }),
                         edges_copy.end());
    } CMR_CATCH_PRINT_THROW("Prepping temp Instance data", err);

    Data::Instance sparse_inst;
    vector<int> elist;
    vector<int> ecap;

    try {

        Graph::get_elist(edges_copy, elist, ecap);
        sparse_inst = Data::Instance("", instance.seed(), ncount, elist, ecap);
    } CMR_CATCH_PRINT_THROW("building sparse Inst", err);

    CCrandstate rstate;
    CCutil_sprand(instance.seed(), &rstate);

    double val = 0;
    int stallcount = std::max(ncount, 250);
    int kicks = std::min(100, (std::max(500, ncount / 2)));

    try { tour.resize(ncount); } CMR_CATCH_PRINT_THROW("allocating cyc", err);

    if (CClinkern_tour(ncount, sparse_inst.ptr(),
                       ecap.size(), &elist[0],
                       stallcount, kicks,
                       const_cast<int *>(&start_tour_nodes[0]),
                       &tour[0], &val,
                       1, 0, 0,
                       (char *) NULL, CC_LK_GEOMETRIC_KICK, &rstate))
        throw runtime_error("CClinkern_tour failed in branch_tour.");

    tour_val = 0.0;
    int fixed_found = 0;

    for (int i = 0; i < ncount; ++i) {
        EndPts e(tour[i], tour[(i + 1) % ncount]);
        tour_val += instance.edgelen(e.end[0], e.end[1]);

        char stat = inds_table(e.end[0], e.end[1]);
        if (stat == 'A') {
            cout << "Avoided edge " << e << " in tour, sparse dummy length "
                 << val << endl;
            found_tour = false;
            return;
        } else if (stat == 'W') {
            ++fixed_found;
        }
    }

    if (fixed_found != want_inds.size()) {
        cout << "Found " << fixed_found << " fixed edges, expected "
             << want_inds.size() << endl;
        found_tour = false;
        return;
    }

    cout << "\tComputed genuine branch tour of length "
         << tour_val << "\n\t"
         << "branch tour: best_tour\t"
         << (tour_val / best_data.min_tour_value) << endl;
}

Sep::Clique::Ptr Executor::compress_tour(const vector<int> &tour)
{
    return tour_cliques.add_tour_clique(tour);
}

vector<int> Executor::expand_tour(const Sep::Clique::Ptr &tour_clique)
{
    return tour_clique->node_list(tour_cliques.ref_tour());
}

void Executor::clamp(const BranchNode &current_node)
{
    int index = core_graph.find_edge_ind(current_node.ends.end[0],
                                         current_node.ends.end[1]);
    if (index == -1)
        throw runtime_error("Tried to clamp edge not in graph.");

    double clamp_bd = static_cast<int>(current_node.direction);
    double clamp_sense = '\0';

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
}

void Executor::unclamp(const BranchNode &current_node)
{
    int index = core_graph.find_edge_ind(current_node.ends.end[0],
                                         current_node.ends.end[1]);
    if (index == -1)
        throw runtime_error("Tried to unclamp edge not in graph.");

    double clamp_bd = static_cast<int>(current_node.direction);
    double clamp_sense;

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
}



}
}
