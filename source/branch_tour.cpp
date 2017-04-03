#include "branch_tour.hpp"
#include "branch_util.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
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

constexpr int IntMax = std::numeric_limits<int>::max();

BranchTourFind::BranchTourFind(const Data::Instance &inst,
                               Data::BestGroup &bestdata,
                               Graph::CoreGraph &coregraph,
                               LP::CoreLP &corelp) try
    : tsp_inst(inst), best_data(bestdata), core_graph(coregraph),
      core_lp(corelp),
      fix_degrees(inst.node_count(), 0),
      tour_edge_tracker(10 * inst.node_count())
{
    int ncount = tsp_inst.node_count();
    int gen_ecount = 0;

    int norm = tsp_inst.ptr()->norm;
    int normbits = (norm & CC_NORM_BITS);

    CCedgegengroup plan;
    CCrandstate rstate;

    CCutil_sprand(tsp_inst.seed(), &rstate);
    CCedgegen_init_edgegengroup(&plan);

    plan.nearest = 4;

    if (norm == CC_EUCLIDEAN || norm == CC_EUCLIDEAN_CEIL) {
        plan.delaunay = 1;
    } else if (normbits == CC_KD_NORM_TYPE || normbits == CC_X_NORM_TYPE) {
        plan.quadnearest = 4;
    } else if (normbits == CC_JUNK_NORM_TYPE) {
        plan.nearest = 12;
    }

    int *elist = NULL;
    if (CCedgegen_edges(&plan, ncount, inst.ptr(), NULL, &gen_ecount,
                        &elist, 1, &rstate))
        throw runtime_error("CCedgegen_edges failed");

    util::c_array_ptr<int> edge_handle(elist);

    extra_edges.reserve(gen_ecount);

    for (int i = 0; i < gen_ecount; ++i) {
        int e0 = elist[2 * i];
        int e1 = elist[(2 * i) + 1];
        if (core_graph.find_edge_ind(e0, e1) == -1)
            extra_edges.emplace_back(e0, e1, tsp_inst.edgelen(e0, e1));
    }

    vector<int> lengths;
    int lengths_size = core_graph.edge_count() + extra_edges.size();

    lengths.reserve(lengths_size);

    for (const Graph::Edge &e : core_graph.get_edges()) {
        lengths.push_back(e.len);
        tour_edge_tracker.add(e.end[0], e.end[1],
                              static_cast<int>(EdgeStats::Core));
    }

    for (const Graph::Edge &e : extra_edges) {
        lengths.push_back(e.len);
        tour_edge_tracker.add(e.end[0], e.end[1],
                              static_cast<int>(EdgeStats::Supply));
    }

    default_length = *std::max_element(lengths.begin(), lengths.end());
    large_length = large_len(core_graph.node_count(), core_graph.get_edges());
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("BranchTourFind constructor failed");
}

/**
 * Computes a branch tour for use as the active tour when cutting and pivoting
 * on \p B. Missing edges will be added to the LP, with tourless mode engaged
 * if no tour can be found.
 * @param[in] B the branch node to be examined next.
 * @param[out] found_tour was a tour found.
 * @param[out] tour if found, these are the tour nodes.
 */
void BranchTourFind::instate_branch_tour(const BranchNode &B,
                                         bool &found_tour, vector<int> &tour)
{
    runtime_error err("Problem in BranchTourFind::instate_branch_tour");

    vector<EndsDir> constraints;

    try { constraints = branch_constraints(B); }
    CMR_CATCH_PRINT_THROW("getting constraints", err);

    bool feas = false;
    double tour_val;

    try {
        compute_tour(constraints, found_tour, feas, tour, tour_val, true);
    } CMR_CATCH_PRINT_THROW("computing tour for estimate", err);

    if (!feas)
        throw runtime_error("Searched for tour on node that should have been "
                            "marked infeasible");

    if (!found_tour) {
        cout << "No tour found, putting core in tourless mode" << endl;
        try { core_lp.tourless_mode(); }
        CMR_CATCH_PRINT_THROW("doing tourless", err);
        return;
    }

    vector<Graph::Edge> missing_edges;
    int ncount = tsp_inst.node_count();
    EdgeStats addstat = EdgeStats::Added;
    if (tour_val < best_data.min_tour_value)
        addstat = EdgeStats::Core;

    try {
        for (int i = 0; i < ncount; ++i) {
            EndPts e(tour[i], tour[(i + 1) % ncount]);
            int e0 = e.end[0];
            int e1 = e.end[1];
            int ind = core_graph.find_edge_ind(e0, e1);
            if (ind == -1) {
                missing_edges.emplace_back(e0, e1, tsp_inst.edgelen(e0, e1));
                tour_edge_tracker.set(e0, e1,
                                      static_cast<int>(addstat));
            }
        }

        int num_new = missing_edges.size();
        if (num_new > 0) {
            if (verbose)
                cout << num_new << " new edges in branch tour" << endl;
            core_lp.add_edges(missing_edges, false);
        }
    } CMR_CATCH_PRINT_THROW("finding/adding missing edges", err);
}

/**
 * Computes a branch tour for use as a pure numeric estimate, without returning
 * the tour itself.
 * @param[in] constraints the constraints the tour must satisfy.
 * @param[out] feas are the constraints free of obvious infeasibilities.
 * @param[out] tour_val the tour length estimate to be set.
 */
void BranchTourFind::estimate_tour(const vector<EndsDir> &constraints,
                                   bool &feas, double &tour_val) try
{
    bool found_tour;
    vector<int> _tour_unused;

    compute_tour(constraints, found_tour, feas, _tour_unused, tour_val,
                 false);
} catch (const exception &e) {
    cerr << e.what() << " computing tour for estimate" << endl;
    throw runtime_error("BranchTourFind::estimate_tour failed");
}

/**
 * To avoid the core edge set becoming extremely large, this function can be
 * called to delete all edges from BranchTourFind#core_graph which are marked
 * in BranchTourFind#tour_edge_tracker as having been added in computing a
 * branch tour. Before doing this, it will mark all edges in the best tour as
 * core graph edges, so as not to accidentally delete edges that helped find
 * an augmenting tour.
 */
void BranchTourFind::prune_edges()
{
    runtime_error err("Problem in BranchTourFind::prune_edges");
    vector<int> delstat;
    int ecount = core_graph.edge_count();

    const vector<int> &best_edges = best_data.best_tour_edges;
    try {
        for (int i = 0; i < best_edges.size(); ++i)
            if (best_edges[i] == 1) {
                const EndPts &e = core_graph.get_edge(i);
                tour_edge_tracker.set(e.end[0], e.end[1],
                                      static_cast<int>(EdgeStats::Core));
            }
    } CMR_CATCH_PRINT_THROW("freezing best tour edges", err);

    try { delstat.resize(0, core_graph.edge_count()); }
    CMR_CATCH_PRINT_THROW("reserving delstat", err)

      for (int i = 0; i < ecount; ++i) {
          const Graph::Edge &e = core_graph.get_edge(i);
          int hval = tour_edge_tracker.get_val(e.end[0], e.end[1]);
          if (hval == EdgeStats::Added)
              delstat[i] = 1;
      }

    try { core_lp.remove_edges(delstat, false); }
    CMR_CATCH_PRINT_THROW("modifying core LP", err);
}

/**
 * @param parent a BranchNode that is currently being split on.
 * @returns a vector of EndsDir objects corresponding to the branchiing
 * constraints for \p parent and all its ancestors.
 * When splitting a problem, we need to compute branch tours for both its
 * children, so we need a list of branch constraints for both of them. This
 * method will return the common branch constraints for the parent node,
 * reserving an extra space in the vector for the current node. Thus if
 * \p parent has child nodes `M` and `N`, to compute a pair of branch tours for
 * `M` and `N` you would call this function on \p parent, and then push/pop
 * the constraints for the individual children.
 */
vector<EndsDir> BranchTourFind::common_constraints(const BranchNode &parent,
                                                   const EndPts &branch_edge)
{
    runtime_error err("Problem in BranchTourFind::common_stats");

    vector<EndsDir> result;
    const BranchNode *iter = &parent;

    try {
        result.reserve(parent.depth + 1);

        while(!iter->is_root()) {
            if (iter->ends == branch_edge) {
                cerr << "Building edge stats for a branch on  "
                     << branch_edge << "\nwith parent "
                     << bnode_brief(parent)
                     << ", encountered " << iter->ends << " as ancestor"
                     << endl;
                throw runtime_error("Trying to split already-branched edge.");
            }

            result.emplace_back(iter->ends, iter->direction);
            iter = iter->parent;
        }
    } CMR_CATCH_PRINT_THROW("building result", err);

    return result;
}

vector<EndsDir> BranchTourFind::branch_constraints(const BranchNode &N)
{
    runtime_error err("Problem in BranchTourFind::branch_constraints");

    vector<EndsDir> result;
    if (!N.is_root())
        try {
            result = common_constraints(*(N.parent), N.ends);
            result.emplace_back(N.ends, N.direction);
        } CMR_CATCH_PRINT_THROW("building result", err);

    return result;
}

/// @returns a Graph::AdjList containing all the edges in \p constraints.
/// Edges fixed to zero will have their AdjObj::val set to 2, or 1 if fixed to
/// one.
Graph::AdjList BranchTourFind::get_fixed_adj(const vector<EndsDir> &constraints)
{
    runtime_error err("Problem in BranchTourFind::get_fixed_adj");

    Graph::AdjList result;
    vector<Graph::Edge> fixed_edges;

    try { fixed_edges.reserve(constraints.size()); }
    CMR_CATCH_PRINT_THROW("reserving fixed edge vec", err);

    for (const EndsDir &e : constraints) {
        const EndPts &ep = e.first;
        int val = 0;
        if (e.second == BranchNode::Dir::Down)
            val = 2;
        else
            val = 1;

        fixed_edges.emplace_back(ep.end[0], ep.end[1], val);
    }

    try { result = Graph::AdjList(tsp_inst.node_count(), fixed_edges); }
    CMR_CATCH_PRINT_THROW("building result", err);

    return result;
}

bool BranchTourFind::tour_compliant(const vector<int> &tour,
                                    const vector<EndsDir> &constraints)
{
    runtime_error err("Problem in BranchTourFind::tour_compliant");

    int ncount = tour.size();
    Graph::AdjList f_adj;

    try { f_adj = get_fixed_adj(constraints); }
    CMR_CATCH_PRINT_THROW("getting fixed adj", err);

    // this loop makes sure no fixed-down edges are in the tour, while marking
    // the fixed-up edges in the tour to check below.
    for (int i = 0; i < ncount; ++i) {
        EndPts e(tour[i], tour[(i + 1) % ncount]);
        int e0 = e.end[0];
        int e1 = e.end[1];

        Graph::AdjObj *e_ptr = f_adj.find_edge(e0, e1);
        if (e_ptr == nullptr) { // edge was not fixed either way
            continue;
        }

        double &e_val = e_ptr->val;

        if (e_val == 2.0) {
            return false; // tour contains down-fixed edge
        } else if (e_val == 1.0) {
            e_val = -1.0; // found a fixed edge in the tour
        }
    }

    // now make sure we hit every fixed-up edge.
    for (const EndsDir &ed : constraints) {
        if (ed.second == BranchNode::Dir::Up) {
            const EndPts &ep = ed.first;
            Graph::AdjObj *e_ptr = f_adj.find_edge(ep.end[0], ep.end[1]);

            if (e_ptr == nullptr) {
                throw err;
            }

            if (e_ptr->val == 1.0) { // then it was not flipped to -1 above
                cout << ep << " is fixed up but not in tour" << endl;
                return false;
            }
        }
    }

    return true;
}

bool BranchTourFind::obvious_infeas(const vector<EndsDir> &constraints)
{
    std::fill(fix_degrees.begin(), fix_degrees.end(), 0);

    for (const EndsDir &ed : constraints)
        for (int pt : ed.first.end)
            if (ed.second == BranchNode::Up) {
                ++(fix_degrees[pt]);
                if (fix_degrees[pt] > 2)
                    return true;
            }

    return false;
}

/**
 * Common method call for trying to find a branch tour based on a set of branch
 * constraints.
 * @param[in] edge_stats the branching constraints.
 * @param[in] start_nodes the starting cycle to use for calling chained LK.
 * @param[out] found_tour was the search able to find a branch tour.
 * @param[out] feas were obvious infeasibilities detected.
 * @param[out] tour the tour nodes vector, if a tour was found.
 * @param[out] tour_val the tour length, if a tour was found.
 * @param[in] bool for_use Is this tour for instatement, or just for a numeric
 * estimate. If true, the BranchTourFind#extra_edges will be added to the
 * sparse instance.
 */
void BranchTourFind::compute_tour(const vector<EndsDir> &edge_stats,
                                  bool &found_tour, bool &feas,
                                  vector<int> &tour, double &tour_val,
                                  bool for_use)
{
    runtime_error err("Problem in BranchTourFind::compute_tour");

    found_tour = true;
    feas = true;
    // if no tour or infeasible tour, output tour vec/length will be voided.
    auto tour_guard = util::make_guard([&found_tour, &feas, &tour, &tour_val]
                                       () -> void
                                       {
                                           if (!found_tour || !feas) {
                                               found_tour = false;
                                               tour.clear();
                                               tour_val = IntMax;
                                           }
                                       });

    try {
        if (obvious_infeas(edge_stats)) {
            feas = false;
            return;
        }
    } CMR_CATCH_PRINT_THROW("checking constraints for infeasibility", err);

    int best_contra_count = 0;
    int active_contra_count = 0;
    const vector<int> &best_tour_edges = best_data.best_tour_edges;
    const vector<double> &active_tour_edges = core_lp.get_active_tour().edges();

    vector<int> want_inds;
    vector<int> avoid_inds;

    try {
        for (const EndsDir &ed : edge_stats) {
            const EndPts &e = ed.first;
            int ind = core_graph.find_edge_ind(e.end[0], e.end[1]);
            if (ind == -1) {
                cerr << "Edge " << e << " not in graph for branch tour" << endl;
                throw err;
            }

            int target_entry = static_cast<int>(ed.second);

            if (target_entry == 0)
                avoid_inds.push_back(ind);
            else
                want_inds.push_back(ind);

            if (best_tour_edges[ind] != target_entry)
                ++best_contra_count;

            if (fabs(active_tour_edges[ind] - target_entry) >= Epsilon::Zero)
                ++active_contra_count;
        }
    } CMR_CATCH_PRINT_THROW("scanning fixed edges", err);

    if (best_contra_count == 0) {
        if (verbose)
            cout << "\tFixed edges affirm best tour, returning it." << endl;

        try { tour = best_data.best_tour_nodes; }
        CMR_CATCH_PRINT_THROW("copying best tour", err);

        tour_val = best_data.min_tour_value;
        return;
    }

    if (active_contra_count == 0) {
        if (verbose)
            cout << "\tFixed edges affirm active tour, returniing it." << endl;
        try { tour = core_lp.get_active_tour().nodes(); }
        CMR_CATCH_PRINT_THROW("copying active tour", err);

        tour_val = core_lp.active_tourlen();
        return;
    }

    const vector<int> &start_tour_nodes =
    ((best_contra_count < active_contra_count) ? best_data.best_tour_nodes :
     core_lp.get_active_tour().nodes());

    vector<Graph::Edge> edges_copy;
    vector<int> elist;
    vector<int> ecap;
    Data::Instance sparse_inst;

    int ncount = tsp_inst.node_count();

    try {
        edges_copy = core_graph.get_edges();
        for (int ind : want_inds)
            edges_copy[ind].len = -large_length;
        for (int ind : avoid_inds)
            edges_copy[ind].len = large_length;

        if (for_use)
            for (const Graph::Edge &e : extra_edges)
                if (tour_edge_tracker.get_val(e.end[0], e.end[1]) == -1)
                    edges_copy.push_back(e);

        Graph::get_elist(edges_copy, elist, ecap);
        sparse_inst = Data::Instance("", tsp_inst.seed(), ncount, elist, ecap,
                                     default_length);

        tour.resize(ncount);
    } CMR_CATCH_PRINT_THROW("prepping data for LK call", err);

    CCrandstate rstate;
    CCutil_sprand(tsp_inst.seed(), &rstate);

    double val = 0;
    int stallcount = std::max(ncount, 250);
    int kicks = std::min(200, (std::max(1000, ncount / 2)));

    if (CClinkern_tour(ncount, sparse_inst.ptr(), ecap.size(), &elist[0],
                       stallcount, kicks,
                       const_cast<int *>(&start_tour_nodes[0]), &tour[0], &val,
                       1, 0, 0, (char *) NULL,
                       CC_LK_CLOSE_KICK, &rstate)) {
        cerr << "CClinkern_tour failed" << endl;
        throw err;
    }

    try {
        found_tour = tour_compliant(tour, edge_stats);
        if (!found_tour)
            return;
    } CMR_CATCH_PRINT_THROW("checking tour compliance", err);

    tour_val = 0.0;

    for (int i = 0; i < ncount; ++i)
        tour_val += tsp_inst.edgelen(tour[i], tour[(i + 1) % ncount]);

    if (verbose) {
        if (for_use)
            cout << "\tComputed genuine branch tour of length "
                 << tour_val << "\n\t"
                 << "branch tour: best_tour\t"
                 << (tour_val / best_data.min_tour_value) << endl;
        else
            cout << "\tComputed branch tour, setting estimate to "
                 << tour_val << "\n\t";
    }

}

}
}
