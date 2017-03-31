#include "branch_tour.hpp"
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

BranchTourFind::BranchTourFind(const Data::Instance &inst,
                               Data::BestGroup &bestdata,
                               Graph::CoreGraph &coregraph) try
    : tsp_inst(inst), best_data(bestdata), core_graph(coregraph),
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

    plan.nearest = 4; //maybe overkill but also could be a helpful safeguard.

    if (norm == CC_EUCLIDEAN || norm == CC_EUCLIDEAN_CEIL) {
        plan.delaunay = 1;
    } else if (normbits == CC_KD_NORM_TYPE || normbits == CC_X_NORM_TYPE) {
        plan.quadnearest = 3;
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

    for (const Graph::Edge &e : core_graph.get_edges())
        lengths.push_back(e.len);

    for (const Graph::Edge &e : extra_edges)
        lengths.push_back(e.len);

    std::nth_element(lengths.begin(), lengths.begin() + (lengths_size / 2),
                     lengths.end());
    default_length = lengths[lengths_size / 2];
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("BranchTourFind constructor failed");
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

}
}
