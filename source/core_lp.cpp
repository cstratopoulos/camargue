#include "core_lp.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <cmath>

using std::cout;
using std::endl;
using std::cerr;

using std::vector;

using std::min;
using std::max;

using std::runtime_error;
using std::logic_error;
using std::exception;

using lpcut_in = CCtsp_lpcut_in;

namespace CMR {
namespace LP {

namespace Eps = CMR::Epsilon;

CoreLP::CoreLP(Data::GraphGroup &graph_data_,
               Data::BestGroup &best_data_) try :
    graph_data(graph_data_), best_data(best_data_),
    ext_cuts(best_data.best_tour_nodes, best_data.perm)
{
    GraphUtils::CoreGraph &core_graph = graph_data.core_graph;
    int ncount = core_graph.node_count();

    char degree_sense = 'E';
    double degree_rhs = 2.0;

    for (int i = 0; i < ncount; ++i)
        new_row(degree_sense, degree_rhs);

    vector<double> col_coeffs{1.0, 1.0};
    double lb = 0.0;
    double ub = 1.0;
    
    for (const CMR::Edge &e : core_graph.get_edges()) {
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

TourBasis::TourBasis(const GraphUtils::CoreGraph &core_graph,
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

double CoreLP::opt()
{
    primal_opt();
    double result = get_objval();
    pivot_back();
    return result;
}

PivType CoreLP::primal_pivot()
{
    runtime_error err("Problem in CoreLP::primal_pivot.");

    double low_limit = best_data.min_tour_value - Eps::Zero;
    GraphUtils::CoreGraph &core_graph = graph_data.core_graph;

    try {
        get_base(tour_base.colstat, tour_base.rowstat);
    } CMR_CATCH_PRINT_THROW("getting base before pivoting", err);


    try {
        nondegen_pivot(low_limit);
        
        get_x(lp_edges);
        
        supp_data.reset(core_graph.node_count(), core_graph.get_edges(),
                        lp_edges, graph_data.island);
        
    } CMR_CATCH_PRINT_THROW("pivoting and setting x", err);

    bool integral = supp_data.integral;
    bool connected = supp_data.connected;
    PivType result;

    if (integral) {
        if (connected) {
            result =  dual_feas() ? PivType::FathomedTour : PivType::Tour;
        } else {
            result =  PivType::Subtour;
        }
    } else {
        result = PivType::Frac;
    }

    if (result == PivType::Tour) {
        try {
            handle_aug();
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

void CoreLP::handle_aug()
{
    double objval = 0;

    best_data.best_tour_nodes = graph_data.island;
    tour_base.best_tour_edges = lp_edges;

    for (int i = 0; i < lp_edges.size(); ++i)
        if (lp_edges[i] < Eps::Zero)
            best_data.best_tour_edges[i] = 0;
        else {
            best_data.best_tour_edges[i] = 1;
            objval += graph_data.core_graph.get_edge(i).len;
        }

    if (objval > best_data.min_tour_value)
        throw runtime_error("Tried to update best tour with worse objval!");

    if (fabs(objval - get_objval()) > Eps::Zero)
        throw runtime_error("Disagreement in new best tour objval with lp.");

    best_data.min_tour_value = objval;

    vector<int> &perm = best_data.perm;
    vector<int> &tour = best_data.best_tour_nodes;
    int ncount = tour.size();

    for (int i = 0; i < ncount; ++i)
        perm[tour[i]] = i;
    
    get_base(tour_base.colstat, tour_base.rowstat);

    vector<double> slacks = row_slacks(ncount, num_rows() - 1);
    
    vector<int> delrows(num_rows(), 0);

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

void CoreLP::add_cuts(Sep::LPcutList &cutq) try
{
    if (cutq.empty())
        return;

    Sep::CutTranslate translator(graph_data);
    vector<int> &perm = best_data.perm;
    vector<int> &tour = best_data.best_tour_nodes;

    for (lpcut_in *cur = cutq.begin(); cur; cur = cur->next) {
        vector<int> rmatind;
        vector<double> rmatval;
        char sense;
        double rhs;

        try {
            translator.get_sparse_row(*cur, perm, rmatind, rmatval, sense,
                                      rhs);
            add_cut(rhs, sense, rmatind, rmatval);
            
        } catch (const exception &e) {
            cerr << e.what() << " adding sparse row.\n";
            throw e;
        }

        try {
            ext_cuts.add_cut(*cur, tour);
        } catch (const exception &e) {
            cerr << e.what() << " adding external rep.\n";
            throw e;
        }
    }
} catch (const exception &e) {
    throw runtime_error("Problem in CoreLP::add_cuts(LPcutList.)");
}

void CoreLP::add_cuts(Sep::CutQueue<Sep::dominoparity> &dpq)
{
    if (dpq.empty())
        return;
    
    runtime_error err("Problem in CoreLP::add_cuts(Sep::dominoparity)");

    Sep::CutTranslate translator(graph_data);
    vector<int> &tour_nodes = best_data.best_tour_nodes;

    for (Sep::CutQueue<Sep::dominoparity>::Itr it = dpq.begin();
         it != dpq.end(); ++it) {
        Sep::dominoparity &dp_cut = *it;
        vector<int> rmatind;
        vector<double> rmatval;
        char sense;
        double rhs;

        try {
            translator.get_sparse_row(dp_cut, tour_nodes, rmatind, rmatval,
                                      sense, rhs);
            add_cut(rhs, sense, rmatind, rmatval);
            ext_cuts.add_cut(dp_cut, rhs, tour_nodes);
        } CMR_CATCH_PRINT_THROW("adding dpcut row", err);
    }
}

void CoreLP::add_edges(vector<EndPts> &add_batch)
{

}

}
}

