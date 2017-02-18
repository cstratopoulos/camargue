#include "core_lp.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>

#include <cmath>

using std::cout;
using std::endl;
using std::cerr;

using std::vector;

using std::abs;
using std::min;
using std::max;

using std::runtime_error;
using std::logic_error;
using std::exception;

using lpcut_in = CCtsp_lpcut_in;


namespace CMR {

using CutType = Sep::HyperGraph::Type;
namespace Eps = CMR::Epsilon;

namespace LP {

CoreLP::CoreLP(Graph::CoreGraph &core_graph_,
               Data::BestGroup &best_data_) try :
    core_graph(core_graph_), best_data(best_data_),
    ext_cuts(best_data.best_tour_nodes, best_data.perm)
{
    int ncount = core_graph.node_count();

    char degree_sense = 'E';
    double degree_rhs = 2.0;

    for (int i = 0; i < ncount; ++i)
        new_row(degree_sense, degree_rhs);

    vector<double> col_coeffs{1.0, 1.0};
    double lb = 0.0;
    double ub = 1.0;

    for (const Graph::Edge &e : core_graph.get_edges()) {
        double objval = e.len;
        vector<int> ends{e.end[0], e.end[1]};

        add_col(objval, ends, col_coeffs, lb, ub);
    }

    tour_base = TourBasis(core_graph, best_data);

    get_row_infeas(tour_base.best_tour_edges,
                   feas_stat, 0, num_rows() - 1);

    copy_start(tour_base.best_tour_edges, tour_base.colstat,
               tour_base.rowstat);

    factor_basis();

    get_x(lp_edges);
    double objval = get_objval();

    if (objval != best_data.min_tour_value) {
        cerr << "Objval after construction: " << objval << "\n";
        throw logic_error("Mismatched obj val.");
    }

    if (lp_edges != tour_base.best_tour_edges)
        throw logic_error("Mismatched lp solution vec.");
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Problem in CoreLP constructor.");
}

TourBasis::TourBasis(const Graph::CoreGraph &core_graph,
                     const Data::BestGroup &best_data) try :
    best_tour_edges(vector<double>(core_graph.edge_count())),
    colstat(vector<int>(core_graph.edge_count(), BStat::AtLower)),
    rowstat(vector<int>(core_graph.node_count(), BStat::AtLower))
{
    const vector<int> &int_tour_edges = best_data.best_tour_edges;

    for(int i = 0; i < int_tour_edges.size(); ++i)
        best_tour_edges[i] = int_tour_edges[i];

    int ncount = core_graph.node_count();
    const vector<int> &tour_nodes = best_data.best_tour_nodes;

    for (int i = 0; i < ncount; ++i) {
        int e0 = tour_nodes[i];
        int e1 = tour_nodes[(i + 1) % ncount];

        int find_ind = core_graph.find_edge_ind(e0, e1);

        if (find_ind == -1) {
            cerr << e0 << ", " << e1 << " not in core graph.\n";
            throw logic_error("Graph does not contain all edges in tour.");
        }

        colstat[find_ind] = BStat::Basic;
    }

    if ((ncount % 2) == 0) {
        int e0 = tour_nodes[ncount - 2];
        int e1 = tour_nodes[ncount - 1];

        int find_ind = core_graph.find_edge_ind(e0, e1);

        if (find_ind == -1) {
            cerr << e0 << ", " << e1 << " not in core graph.\n";
            throw logic_error("Graph does not contain all edges in tour.");
        }

        colstat[find_ind] = BStat::AtUpper;

        e0 = tour_nodes[0];
        e1 = tour_nodes[ncount - 2];

        if (find_ind == -1) {
            cerr << e0 << ", " << e1 << " not in core graph.\n";
            throw logic_error("Graph does not contain all edges in tour.");
        }

        colstat[find_ind] = BStat::Basic;
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("TourBasis constructor failed.");
}

PivType CoreLP::primal_pivot()
{
    runtime_error err("Problem in CoreLP::primal_pivot");

    double low_limit = best_data.min_tour_value - Eps::Zero;

    try {
        get_base(tour_base.colstat, tour_base.rowstat);
    } CMR_CATCH_PRINT_THROW("getting base before pivoting", err);

    vector<int> dfs_island;

    try {
        nondegen_pivot(low_limit);
        get_x(lp_edges);
        supp_data = Data::SupportGroup(core_graph.get_edges(),
                                       lp_edges, dfs_island,
                                       core_graph.node_count());
    } CMR_CATCH_PRINT_THROW("pivoting and setting x", err);

    ++num_nd_pivots;
    sum_it_count += it_count();

    bool integral = supp_data.integral;
    bool connected = supp_data.connected;
    PivType result = PivType::Frac;

    if (integral) {
        if (connected) {
            result =  dual_feas() ? PivType::FathomedTour : PivType::Tour;
        } else {
            result =  PivType::Subtour;
        }
    } else {
        if (dual_feas()) {
            if (get_objval() >= best_data.min_tour_value - 1.0 + Eps::Zero)
                result = PivType::FathomedTour;
        } else
            result = PivType::Frac;
    }

    if (result == PivType::Tour) {
        try {
            handle_aug_pivot(dfs_island);
        } CMR_CATCH_PRINT_THROW("handling augmentation", err);
    }

    return result;
}

void CoreLP::pivot_back()
{
    try {
        copy_start(tour_base.best_tour_edges, tour_base.colstat,
                   tour_base.rowstat);
        factor_basis();
    } catch (const exception &e) {
        cerr << e.what() << "\n";
        throw runtime_error("Problem in CoreLP::pivot_back.");
    }
}

void CoreLP::handle_aug_pivot(const std::vector<int> &tour_nodes)
{
    runtime_error err("Problem in CoreLP::handle_aug_pivot");

    best_data.best_tour_nodes = tour_nodes;
    tour_base.best_tour_edges = lp_edges;

    for (int i = 0; i < lp_edges.size(); ++i)
        if (lp_edges[i] < Eps::Zero)
            best_data.best_tour_edges[i] = 0;
        else
            best_data.best_tour_edges[i] = 1;

    try {
        update_best_data();
        get_base(tour_base.colstat, tour_base.rowstat);
    } CMR_CATCH_PRINT_THROW("instating tour vec", err);

    try { prune_slacks(); } CMR_CATCH_PRINT_THROW("pruning slacks", err);
}

void CoreLP::set_best_tour(const std::vector<int> &tour_nodes)
{
    runtime_error err("Problem in set_best_tour");
    best_data.best_tour_nodes = tour_nodes;

    int ncount = tour_nodes.size();
    vector<int> &tour_edges = best_data.best_tour_edges;
    vector<double> &d_tour_edges = tour_base.best_tour_edges;

    for (int i = 0; i < tour_edges.size(); ++i) {
        tour_edges[i] = 0;
        d_tour_edges[i] = 0.0;
    }

    for (int i = 0; i < ncount; ++i) {
        EndPts e(tour_nodes[i], tour_nodes[(i + 1) % ncount]);
        int ind = core_graph.find_edge_ind(e.end[0], e.end[1]);
        if (ind == -1) {
            cerr << "Edge " << e.end[0] << ", " << e.end[1]
                 << " still not in graph\n";
            throw err;
        }
        tour_edges[ind] = 1;
        d_tour_edges[ind] = 1.0;
    }

    try {
        update_best_data();
        copy_start(d_tour_edges);
        factor_basis();
        get_base(tour_base.colstat, tour_base.rowstat);
    } CMR_CATCH_PRINT_THROW("instating tour vec", err);

    try { prune_slacks(); } CMR_CATCH_PRINT_THROW("pruning slacks", err);
}

void CoreLP::update_best_data()
{
    double edge_objval = 0.0;

    const vector<int> &tour_edges = best_data.best_tour_edges;

    for (int i = 0; i < tour_edges.size(); ++i)
        if (tour_edges[i] == 1)
            edge_objval += core_graph.get_edge(i).len;

    if (edge_objval > best_data.min_tour_value) {
        cerr << "Edge objval is " << edge_objval << "\n";
        throw runtime_error("Tried to update best tour with worse objval!");
    }

    best_data.min_tour_value = edge_objval;

    vector<int> &perm = best_data.perm;
    vector<int> &tour = best_data.best_tour_nodes;
    int ncount = tour.size();

    for (int i = 0; i < ncount; ++i)
        perm[tour[i]] = i;
}

void CoreLP::prune_slacks()
{
    get_x(lp_edges);

    for (int i = 0; i < lp_edges.size(); ++i)
        if (fabs(lp_edges[i] - tour_base.best_tour_edges[i]) >= Eps::Zero)
            throw runtime_error("Tried to prune slacks with non-tour vec");

    int ncount = core_graph.node_count();

    vector<double> slacks = row_slacks(ncount, num_rows() - 1);

    int orig_numrows = num_rows();
    vector<int> delrows(orig_numrows, 0);

    int rownum = ncount;

    for (double slack : slacks) {
        if (slack)
            delrows[rownum] = 1;
        ++rownum;
    }

    del_set_rows(delrows);
    ext_cuts.del_cuts(delrows);
    factor_basis();
}

void CoreLP::rebuild_basis()
{
    vector<double> &tour = tour_base.best_tour_edges;

    copy_start(tour);
    factor_basis();

    vector<double> test_x = lp_vec();

    for (int i = 0; i < test_x.size(); ++i)
        if (std::abs(test_x[i] - tour[i]) >= Eps::Zero)
            throw runtime_error("tour not instated after basis rebuild.");

    get_base(tour_base.colstat, tour_base.rowstat);

    vector<double> feas_stat;

    get_row_infeas(tour, feas_stat, 0, num_rows() - 1);

    int ncount = core_graph.node_count();

    for (int i = 0; i < feas_stat.size(); ++i) {
        if (std::abs(feas_stat[i] >= Eps::Zero)) {
            if (i < ncount)
                cout << "Found nonzero infeas on degree eqn\n";
            else {
                cout << "Found nonzero infeas on cut, type: "
                     << ext_cuts.get_cut(i).cut_type() << "\n";
            }

            throw runtime_error("tour is now infeasible.");
        }
    }
}

void CoreLP::add_cuts(Sep::LPcutList &cutq)
{
    if (cutq.empty())
        return;

    runtime_error err("Problem in CoreLP::add_cuts LPcutList");

    vector<int> &perm = best_data.perm;
    vector<int> &tour = best_data.best_tour_nodes;
    try {
        for (const lpcut_in *cur = cutq.begin(); cur; cur = cur->next) {
            SparseRow R = Sep::get_row(*cur, perm, core_graph);
            add_cut(R);
            ext_cuts.add_cut(*cur, tour);
        }
    } CMR_CATCH_PRINT_THROW("processing/adding cut", err);

    cutq.clear();
}

void CoreLP::add_cuts(Sep::CutQueue<Sep::dominoparity> &dpq)
{
    runtime_error err("Problem in CoreLP::add_cuts(Sep::dominoparity)");

    vector<int> &tour_nodes = best_data.best_tour_nodes;

    try {
        while(!dpq.empty()) {
            const Sep::dominoparity &dp_cut = dpq.peek_front();
            SparseRow R = Sep::get_row(dp_cut, tour_nodes, core_graph);
            add_cut(R);
            ext_cuts.add_cut(dp_cut, R.rhs, tour_nodes);
            dpq.pop_front();
        }
    } CMR_CATCH_PRINT_THROW("processing/adding cut", err);
}

void CoreLP::add_cuts(Sep::CutQueue<SparseRow> &gmi_q)
{
    runtime_error err("Problem in CoreLP::add_cuts(Sep::SparseRow)");

    try {
        while (!gmi_q.empty()) {
            add_cut(gmi_q.peek_front());
            ext_cuts.add_cut();
            gmi_q.pop_front();
        }
    } CMR_CATCH_PRINT_THROW("adding sparse cut row", err);
}

void CoreLP::add_cuts(Sep::CutQueue<Sep::ex_blossom> &ex2m_q)
{
    runtime_error err("Problem in CoreLP::add_cuts(Sep::ex_blossom)");

    int ncount = core_graph.node_count();

    try {
        while (!ex2m_q.empty()) {
            const Sep::ex_blossom &B = ex2m_q.peek_front();
            const vector<int> &handle_nodes = B.handle;
            const vector<Graph::Edge> &edges = core_graph.get_edges();

            vector<int> handle_delta  = Graph::delta_inds(handle_nodes, edges,
                                                          ncount);
            vector<int> tooth_inds = Sep::teeth_inds(B,
                                                     best_data.best_tour_edges,
                                                     lp_edges, edges, ncount,
                                                     handle_delta);
            vector<vector<int>> tooth_edges;
            tooth_edges.reserve(tooth_inds.size());

            for (int ind : tooth_inds) {
                const Graph::Edge &e = edges[ind];
                tooth_edges.emplace_back(vector<int>{e.end[0], e.end[1]});
            }

            SparseRow R = Sep::get_row(handle_delta, tooth_edges, core_graph);
            add_cut(R);
            ext_cuts.add_cut(handle_nodes, tooth_edges);
            ex2m_q.pop_front();
        }
    } CMR_CATCH_PRINT_THROW("processing/adding cuts", err);
}

void CoreLP::add_cuts(Sep::CutQueue<Sep::HyperGraph> &pool_q)
{
    runtime_error err("Problem in CoreLP::add_cuts(Sep::HyperGraph)");

    try {
        while (!pool_q.empty()) {
            Sep::HyperGraph &H = pool_q.peek_front();
            SparseRow R;

            R.sense = H.get_sense();
            R.rhs = H.get_rhs();
            H.get_coeffs(core_graph.get_edges(), R.rmatind, R.rmatval);
            add_cut(R);
            ext_cuts.add_cut(H);
            pool_q.pop_front();
        }
    } CMR_CATCH_PRINT_THROW("processing/adding cuts", err);
}

void CoreLP::add_edges(const vector<Graph::Edge> &batch)
{
    runtime_error err("Problem in CoreLP::add_edges");

    int old_ecount = core_graph.edge_count();
    int new_ecount = old_ecount + batch.size();

    try {
        best_data.best_tour_edges.resize(new_ecount, 0);
        for (const Graph::Edge &e : batch)
            core_graph.add_edge(e);
    } CMR_CATCH_PRINT_THROW("adding edges to core graph/best group", err);

    double lb = 0.0;
    double ub = 1.0;
    vector<int> cmatind;
    vector<double> cmatval;

    try {
        for (const Graph::Edge &e : batch) {
            double objval = e.len;

            ext_cuts.get_col(e.end[0], e.end[1], cmatind, cmatval);
            add_col(objval, cmatind, cmatval, lb, ub);
        }
    } CMR_CATCH_PRINT_THROW("adding edges to core lp", err);

    try {
        tour_base.best_tour_edges.resize(new_ecount, 0.0);
        lp_edges.resize(new_ecount, 0.0);
    } CMR_CATCH_PRINT_THROW("resizing for edges", err)
}

void CoreLP::purge_gmi()
{
    vector<int> delrows(num_rows(), 0);
    int ncount = core_graph.node_count();

    int delcount = 0;
    int i = 0;
    for (const Sep::HyperGraph &H : ext_cuts.get_cuts()) {
        if (H.cut_type() == CutType::Non) {
            delrows[ncount + i] = 1;
            ++delcount;
        }
        ++i;
    }

    if (delcount > 0) {
        del_set_rows(delrows);
        ext_cuts.del_cuts(delrows);
        factor_basis();
    }
}

}
}
