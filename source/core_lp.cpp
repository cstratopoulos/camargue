#include "core_lp.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <stdexcept>

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
    active_tour(core_graph_, best_data),
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

    instate_active();

    get_x(lp_edges);

    if (lp_edges != active_tour.edges())
        throw logic_error("Mismatched lp solution vec.");
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Problem in CoreLP constructor.");
}

/**
 * This method is responsible for computing and reporting information about
 * a non-degenerate primal pivot, invoking methods to update bounds, prune cuts,
 * etc., if necessary.
 * @returns a PivType corresponding to the vector obtained from the
 * non-degenerate primal pivot.
 */
PivType CoreLP::primal_pivot()
{
    runtime_error err("Problem in CoreLP::primal_pivot");

    int ncount = core_graph.node_count();

    vector<int> dfs_island;
    Basis bas;

    try {
        bas = basis_obj();
        nondegen_pivot(active_tourlen());
        //cb_nondegen_pivot(active_tourlen(), bas, 1);
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
    bool dfeas = dual_feas();
    PivType result = PivType::Frac;

    if (integral) {
        if (connected) {
            if (dfeas) {
                result = PivType::FathomedTour;
                genuine_opt_tour = true;
            } else {
                result = PivType::Tour;
            }
        } else {
            result =  PivType::Subtour;
        }
    } else {
        if (dfeas) {
            if (get_objval() >= global_ub() - 0.9)
                result = PivType::FathomedTour;
        } else
            result = PivType::Frac;
    }

    if(is_tour_piv(result))
        ext_cuts.reset_ages();
    else if (num_rows() > core_graph.node_count())
        try {
            ext_cuts.piv_age_cuts(pi(core_graph.node_count(), num_rows() - 1));
        } CMR_CATCH_PRINT_THROW("upating pivot ages", err);

    if (result == PivType::Tour) {
        try {
            handle_aug_pivot(std::move(dfs_island), basis_obj());
        } CMR_CATCH_PRINT_THROW("handling augmentation", err);
    } else if (genuine_opt_tour) {
        try {
            active_tour.set_basis(basis_obj());
        } CMR_CATCH_PRINT_THROW("setting opt tour base", err);
    } else {
        if (!bas.empty()) {
            active_tour.set_basis(std::move(bas));
        }
    }

    return result;
}

/**
 * This method instates active_tour.
 * @param prune_slacks if true, this function will delete all cuts whose slack
 * variables are in the basis for active_tour.
 */
void CoreLP::pivot_back(bool prune_slacks)
{
    string error_string =
    "Problem in CoreLP::pivot_back, prune " + std::to_string(prune_slacks);

    runtime_error err(error_string);

    int delct = 0;
    int numrows = num_rows();
    int rowdiff = numrows - prev_numrows;

    const vector<int> &cut_stats = active_tour.base().rowstat;
    vector<int> delset;

    if (prune_slacks) {
        try { delset = vector<int>(numrows, 0); }
        CMR_CATCH_PRINT_THROW("allocating delset", err);

        for (int i = core_graph.node_count(); i < prev_numrows; ++i) {
            const Sep::HyperGraph &H = ext_cuts.get_cut(i);
            if (H.fresh_cut())
                continue;

            if (H.tour_age() >= CutAge::TourOld &&
                H.piv_age() >= CutAge::PivOld) {
                delset[i] = 3;
                ++delct;
            }
        }

        for (int i = prev_numrows; i < numrows; ++i)
            if (cut_stats[i] == 1) {
                delset[i] = 1;
                ++delct;
            }
    }

    if (delct > 0 && delct != rowdiff) {
        try {
            ext_cuts.del_cuts(delset);
            del_set_rows(delset);
            prev_numrows = num_rows();
            reset_instate_active();
        } CMR_CATCH_PRINT_THROW("deleting cuts/instating tour", err);
    } else {
        try {
            instate_active();
        } CMR_CATCH_PRINT_THROW("instating tour", err);
    }

    if (num_rows() > core_graph.node_count())
        try {
            ext_cuts.tour_age_cuts(pi(core_graph.node_count(),
                                      num_rows() - 1));
        } CMR_CATCH_PRINT_THROW("updating tour ages", err);
}

void CoreLP::instate_active()
{
    runtime_error err("CoreLP::instate_active failed");

    try {
        active_tour.instate(*this);
    } catch (const exception &e) {
        cerr << e.what() << " instating active tour" << endl;

        bool active_feas;
        bool best_feas;

        try {
            active_feas = check_feas(active_tour.edges());
            best_feas = check_feas(best_data.best_tour_edges);
        } catch (const exception &e) {
            cerr << e.what() << " getting feas stats" << endl;
            throw err;
        }

        cout << "Active tour vector feasible: " << active_feas << "\n"
             << "Best tour vector feasible: " << best_feas << endl;
        throw err;
    }
}

void CoreLP::reset_instate_active()
{
    runtime_error err("CoreLP::reset_instate_active failed");

    try {
        active_tour.reset_instate(*this);
    } catch (const exception &e) {
        cerr << e.what() << " reset instating active tour" << endl;

        bool active_feas;
        bool best_feas;

        try {
            active_feas = check_feas(active_tour.edges());
            best_feas = check_feas(best_data.best_tour_edges);
        } catch (const exception &e) {
            cerr << e.what() << " getting feas stats" << endl;
            throw err;
        }

        cout << "Active tour vector feasible: " << active_feas << "\n"
             << "Best tour vector feasible: " << best_feas << endl;
        throw err;
    }
}

/**
 * @param tour_nodes the sequence of tour nodes giving a tour which augments
 * active_tour. Should be obtained from the depth-first search that was used
 * to establish connectivity of supp_data.
 * @param aug_base the basis for the augmenting pivot.
 */
void CoreLP::handle_aug_pivot(vector<int> tour_nodes, Basis aug_base)
{
    runtime_error err("Problem in CoreLP::handle_aug_pivot");

    try {
        active_tour = ActiveTour(std::move(tour_nodes),
                                 lp_vec(), std::move(aug_base),
                                 get_objval(), core_graph);
    } catch (const exception &e) {
        cerr << e.what() << " constructing active tour from pivot" << endl;

        try {
            vector<double> x = lp_vec();
            bool augpiv_feas = check_feas(x);
            bool best_feas = check_feas(best_data.best_tour_edges);
            cout << "Augmenting lp vec feasible: " << augpiv_feas << "\n"
                 << "Best tour vector feasible: " << best_feas << endl;
        } catch (const exception &e) {
            cerr << e.what() << " trying to check feasibility" << endl;
            throw err;
        }
        throw err;
    }

    try { prune_slacks(); } CMR_CATCH_PRINT_THROW("pruning slacks", err);
}

/**
 * For setting active_tour by a sequence of tour nodes, presumably obtained
 * from Solver::frac_recover. This function will generate an edge vector
 * representation and a basis, instating the tour.
 */
void CoreLP::set_active_tour(std::vector<int> tour_nodes)
{
    runtime_error err("Problem in CoreLP::set_active_tour");
    vector<int> tnodes_copy;

    try {
        tnodes_copy = tour_nodes;
        active_tour = ActiveTour(std::move(tour_nodes), *this,
                                 core_graph);
    } catch (const exception &e) {
        cerr << e.what() << " constructing active tour from nodes" << endl;

        try {
            vector<int> set_edges;
            double _tval;

            core_graph.tour_edge_vec(tnodes_copy, set_edges, _tval);

            bool best_feas = check_feas(best_data.best_tour_edges);
            bool set_feas = check_feas(set_edges);
            cout << "Would-be tour vector feasible: " << set_feas << "\n"
                 << "Best tour vector feasible: " << best_feas << endl;
        } catch (const exception &e) {
            cerr << e.what() << " trying to check feasibility" << endl;
            throw err;
        }
        throw err;
    }

    try { prune_slacks(); } CMR_CATCH_PRINT_THROW("pruning slacks", err);
}


void CoreLP::prune_slacks()
{
    if (num_rows() == core_graph.node_count())
        return;

    get_x(lp_edges);

    for (int i = 0; i < lp_edges.size(); ++i)
        if (fabs(lp_edges[i] - active_tour.edges()[i]) >= Eps::Zero)
            throw runtime_error("Tried to prune slacks with non-tour vec");

    int ncount = core_graph.node_count();

    vector<double> slacks = row_slacks(ncount, num_rows() - 1);

    int orig_numrows = num_rows();
    vector<int> delrows(orig_numrows, 0);

    int rownum = ncount;
    int delcount = 0;

    for (double slack : slacks) {
        if (slack) {
            delrows[rownum] = 3;
            ++delcount;
        }
        ++rownum;
    }

    if (verbose)
        cout << "\t" << delcount << " slack rows can be pruned" << endl;

    ext_cuts.del_cuts(delrows);
    del_set_rows(delrows);
    prev_numrows = num_rows();

    reset_instate_active();
}

void CoreLP::add_cuts(Sep::LPcutList &cutq)
{
    if (cutq.empty())
        return;
    prev_numrows = num_rows();

    runtime_error err("Problem in CoreLP::add_cuts LPcutList");

    const vector<int> &perm = active_tour.tour_perm();
    const vector<int> &tour = active_tour.nodes();
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

    const vector<int> &tour_nodes = active_tour.nodes();

    try {
        while(!dpq.empty()) {
            const Sep::dominoparity &dp_cut = dpq.peek_front();
            SparseRow R = Sep::get_row(dp_cut, tour_nodes, core_graph);

            // This should be investigated later but sometimes extremely
            // dense cuts are returned, so this is a bad hacky workaround.
            if (R.rmatind.size() < core_graph.edge_count() / 4) {
                add_cut(R);
                ext_cuts.add_cut(dp_cut, R.rhs, tour_nodes);
            }
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

    const vector<double> &tour_edges = active_tour.edges();

    try {
        while (!ex2m_q.empty()) {
            const Sep::ex_blossom &B = ex2m_q.peek_front();
            const vector<int> &handle_nodes = B.handle;
            const vector<Graph::Edge> &edges = core_graph.get_edges();

            vector<int> handle_delta  = Graph::delta_inds(handle_nodes, edges,
                                                          ncount);
            vector<int> tooth_inds = Sep::teeth_inds(B,
                                                     tour_edges, lp_edges,
                                                     edges, ncount,
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

/**
 * @param reinstate if true, method will attempt to reinstate the active_tour.
 */
void CoreLP::add_edges(const vector<Graph::Edge> &batch, bool reinstate)
{
    runtime_error err("Problem in CoreLP::add_edges");

    int old_ecount = core_graph.edge_count();
    int new_ecount = old_ecount + batch.size();

    try {
        purge_gmi(false);
    } CMR_CATCH_PRINT_THROW("getting rid of gmi", err);

    try {
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
        lp_edges.resize(new_ecount, 0.0);
        best_data.best_tour_edges.resize(new_ecount, 0);
    } CMR_CATCH_PRINT_THROW("adding edges to core lp/resizing", err);

    if (reinstate) {
        try { reset_instate_active(); }
        CMR_CATCH_PRINT_THROW("resetting active tour", err);
    }
}

/**
 * @param edge_delstat a vector of length equal to `core_graph.edge_count()`,
 * with `delstat[i]` equal to one if the corresponding edge is to be removed
 * from the core_graph and the LP.
 */
void CoreLP::remove_edges(vector<int> edge_delstat)
{
    runtime_error err("Problem in CoreLP::add_edges");

    try {
        purge_gmi(false);
    } CMR_CATCH_PRINT_THROW("getting rid of gmi", err);

    int ecount = core_graph.edge_count();
    if (edge_delstat.size() != ecount)
        throw runtime_error("Size mismatch in remove_edges");

    vector<Graph::Edge> &graph_edges = core_graph.get_edges();
    vector<int> &delstat = edge_delstat;

    for (int i = 0; i < ecount; ++i)
        if (delstat[i] == 1)
            graph_edges[i].removable = true;
        else
            graph_edges[i].removable = false;

    try {
        core_graph.remove_edges();
    } CMR_CATCH_PRINT_THROW("removing edges from coregraph", err);

    lp_edges.resize(core_graph.edge_count());

    vector<int> &tour_edges = best_data.best_tour_edges;
    for (int i = 0; i < delstat.size(); ++i)
        if (delstat[i] == 1)
            tour_edges[i] = -1;

    tour_edges.erase(std::remove(tour_edges.begin(), tour_edges.end(), -1),
                     tour_edges.end());

    try {
        del_set_cols(delstat);
        active_tour = ActiveTour(best_data.best_tour_nodes, *this,
                                 core_graph);
    } CMR_CATCH_PRINT_THROW("deleting cols from LP/reinstating", err);


}

/**
 * @param instate if true, try to instate the tour after cut deletion.
 */
void CoreLP::purge_gmi(bool instate)
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
        ext_cuts.del_cuts(delrows);
        del_set_rows(delrows);
        prev_numrows = num_rows();

        if (instate)
            reset_instate_active();
    }
}

/**
 * @tparam numtype the number type for the vector to check.
 * @param x_vec the vector to examine.
 * @returns true if \p x_vec satisfies all constraints and column bounds in the
 * LP, false otherwise.
 * @remarks Prints EXTREMELY verbose information if a row infeasibility is
 * found.
 */
template<typename numtype>
bool CoreLP::check_feas(const std::vector<numtype> &x_vec)
{
    using DoublePair = std::pair<double, double>;

    runtime_error err("Problem in CoreLP::check_feas");
    bool result = true;

    int numrows = num_rows();
    int numcols = num_cols();

    vector<double> dbl_x;
    vector<double> row_feas;
    vector<double> col_feas;

    cout << "Reporting if solution is feasible, " << numrows << " rows and "
         << numcols << " cols in the LP\n" << ext_cuts.get_cuts().size()
         << " external cuts" << endl;

    try {
        dbl_x = vector<double>(x_vec.begin(), x_vec.end());
        get_row_infeas(dbl_x, row_feas, 0, numrows - 1);
        get_col_infeas(dbl_x, col_feas, 0, numcols - 1);
    } CMR_CATCH_PRINT_THROW("making solver queries", err);

    int ncount = core_graph.node_count();

    for (int i = 0; i < numrows; ++i) {
        double rowfeas = row_feas[i];
        if (rowfeas != 0.0) {
            result = false;
            cout << "\nRow " << i << " has infeas of "  << rowfeas << " on ";
            if (i < ncount)
                cout << "degree equation" << endl;
            else {
                const Sep::HyperGraph &H = ext_cuts.get_cut(i);
                cout << H.cut_type() << "\n\tsense "
                     << H.get_sense() << ", rhs " << H.get_rhs()
                     << ", clq/tooth count " << H.get_cliques().size() << "/"
                     << H.get_teeth().size() << endl;
                LP::SparseRow rel_row;
                LP::SparseRow hg_row;

                try {
                    rel_row = get_row(i);
                    H.get_coeffs(core_graph.get_edges(), hg_row.rmatind,
                                 hg_row.rmatval);
                } CMR_CATCH_PRINT_THROW("getting SparseRows", err);
                cout << "HG row act: " << Sep::get_activity(dbl_x, hg_row)
                     << ", rel row act: " << Sep::get_activity(dbl_x, rel_row)
                     << endl;
                cout << "HG row size " << hg_row.rmatind.size() << ", rel "
                     << rel_row.rmatind.size() << endl;
                cout << "Rel row sense " << rel_row.sense << ", rhs "
                     << rel_row.rhs << endl;

                vector<int> coeff_compare(numcols, 1);
                vector<DoublePair> hg_rel_coeffs(numcols,
                                                 DoublePair(-1.0, -1.0));

                for (int i = 0; i < hg_row.rmatind.size(); ++i) {
                    coeff_compare[hg_row.rmatind[i]] *= 2;
                    hg_rel_coeffs[hg_row.rmatind[i]].first = hg_row.rmatval[i];
                }

                for (int i = 0; i < rel_row.rmatind.size(); ++i) {
                    coeff_compare[rel_row.rmatind[i]] *= 3;
                    hg_rel_coeffs[rel_row.rmatind[i]].second =
                    rel_row.rmatval[i];
                }

                for (int i = 0; i < numcols; ++i) {
                    double entry = coeff_compare[i];
                    if (entry == 2)
                        cout << "Coeff " << i << " only in HG with val "
                             << static_cast<int>(hg_rel_coeffs[i].first)
                             << "\n";
                    else if (entry == 3)
                        cout << "Coeff " << i << " only in rel with val "
                             << static_cast<int>(hg_rel_coeffs[i].second)
                             << "\n";
                }


            }
        }
    }

    vector<double> lbs;
    vector<double> ubs;

    try {
        lbs = lower_bds(0, numcols - 1);
        ubs = upper_bds(0, numcols - 1);
    } CMR_CATCH_PRINT_THROW("getting bounds", err);

    for (int i = 0; i < numcols; ++i) {
        double colfeas = col_feas[i];
        if (colfeas != 0.0) {
            result = false;
            Graph::Edge e = core_graph.get_edge(i);
            cout << "Found infeas of " << colfeas << " on edge "
                 << e << " with bds " << static_cast<int>(lbs[i])
                 << " <= x_" << i
                 << " <= " << static_cast<int>(ubs[i]) << "\n"
                 << "Endpoint degrees "
                 << core_graph.get_adj().nodelist[e.end[0]].degree()
                 << " / "
                 << core_graph.get_adj().nodelist[e.end[1]].degree()
                 << endl;
        }
    }

    cout << "Report completed with result " << result << endl;
    return result;
}

}
}
