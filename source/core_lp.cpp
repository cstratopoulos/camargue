#include "core_lp.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <string>

#include <cmath>

using std::cout;
using std::endl;
using std::cerr;
using std::string;

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
    ext_cuts(best_data.best_tour_nodes, best_data.perm),
    prev_numrows(core_graph_.node_count())
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

    copy_start(tour_base.best_tour_edges, tour_base.base.colstat,
               tour_base.base.rowstat);

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
    best_tour_edges(vector<double>(core_graph.edge_count()))
{
    const vector<int> &int_tour_edges = best_data.best_tour_edges;
    vector<int> &colstat = base.colstat;
    vector<int> &rowstat = base.rowstat;

    colstat = vector<int>(core_graph.edge_count(), BStat::AtLower);
    rowstat = vector<int>(core_graph.node_count(), BStat::AtLower);

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

    int ncount = core_graph.node_count();

    vector<int> dfs_island;
    Basis bas;

    try {
        //nondegen_pivot(best_data.min_tour_value);
        cb_nondegen_pivot(best_data.min_tour_value, bas, 1);
        get_x(lp_edges);
        supp_data = Data::SupportGroup(core_graph.get_edges(),
                                       lp_edges, dfs_island,
                                       ncount);
    } CMR_CATCH_PRINT_THROW("pivoting and setting x", err);

    ++num_nd_pivots;
    sum_it_count += it_count();

    bool integral = supp_data.integral;
    bool connected = supp_data.connected;
    bool genuine_opt_tour = false;
    PivType result = PivType::Frac;

    if (integral) {
        if (connected) {
            if (dual_feas()) {
                result = PivType::FathomedTour;
                genuine_opt_tour = true;
            } else {
                result = PivType::Tour;
            }
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
    } else if (genuine_opt_tour) {
        try {
            tour_base.base = basis_obj();
        } CMR_CATCH_PRINT_THROW("setting opt tour base", err);
    } else {
        if (!bas.empty()) {
            tour_base.base = std::move(bas);
        }
    }


    // try {
    //     if (is_tour_piv(result))
    //         cut_mon = LP::CutMonitor(num_rows() - ncount);
    //     else {
    //         if (num_rows() > ncount) {
    //             get_pi(pi_vals, ncount, num_rows() - 1);
    //             cut_mon.update_pivs(pi_vals);
    //         }
    //     }
    // } CMR_CATCH_PRINT_THROW("updating CutMonitor", err);

    return result;
}

/**
 * @post the tour in best_data is instated as the current LP solution with
 * an associated row and column basis. Calls to lp_vec, get_objval, etc. are
 * valid and return the best tour.
 * @param prune_slacks if true, this function will delete all cuts whose slack
 * variables are in the current tour basis.
 */
void CoreLP::pivot_back(bool prune_slacks)
{
    string error_string =
    "Problem in CoreLP::pivot_back, prune " + std::to_string(prune_slacks);

    runtime_error err(error_string);

    int delct = 0;
    int numrows = num_rows();
    int rowdiff = numrows - prev_numrows;

    vector<int> &cut_stats = tour_base.base.rowstat;
    vector<int> delset(numrows, 0);

    if (prune_slacks)
        for (int i = prev_numrows; i < numrows; ++i)
            if (cut_stats[i] == 1) {
                delset[i] = 1;
                ++delct;
            }

    if (delct > 0 && delct != rowdiff) {
        try {
            del_set_rows(delset);
            ext_cuts.del_cuts(delset, false);
            copy_start(tour_base.best_tour_edges);
            factor_basis();
            // cout << "Deleted " << delct << " / "
            //      << (numrows - prev_numrows) << " basic slack cuts" << endl;
        } CMR_CATCH_PRINT_THROW("deleting cuts/factoring tour basis", err);
    } else {
        try {
            copy_start(tour_base.best_tour_edges);
            // copy_start(tour_base.best_tour_edges, tour_base.colstat,
            //            tour_base.rowstat);
            factor_basis();
        } CMR_CATCH_PRINT_THROW("copying/factoring basis", err);
    }

    double objval = get_objval();
    if (fabs(objval - best_tourlen()) >= Eps::Zero) {
        cerr << "Objval disagreement: " << objval << " vs tour length "
             << best_tourlen() << endl;
        throw err;
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
        tour_base.base = basis_obj();
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

    std::fill(tour_edges.begin(), tour_edges.end(), 0);
    std::fill(d_tour_edges.begin(), d_tour_edges.end(), 0.0);

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
        tour_base.base = basis_obj();
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
    ext_cuts.del_cuts(delrows, true);
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

    tour_base.base = basis_obj();

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
    prev_numrows = num_rows();

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
    prev_numrows = num_rows();

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
    prev_numrows = num_rows();

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
    prev_numrows = num_rows();

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
    prev_numrows = num_rows();

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
        ext_cuts.del_cuts(delrows, false);
        factor_basis();
    }
}

}
}
