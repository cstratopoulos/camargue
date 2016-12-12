#include "core_lp.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <vector>
#include <iostream>
#include <stdexcept>

using std::cout;
using std::cerr;

using std::vector;

using std::min;
using std::max;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {
namespace LP {

CoreLP::CoreLP(CMR::Data::GraphGroup &graph_data_,
               CMR::Data::BestGroup &best_data_) try :
    graph_data(graph_data_), best_data(best_data_)
{
    CMR::Graph &graph = graph_data.m_graph;
    int ncount = graph.node_count;

    char degree_sense = 'E';
    double degree_rhs = 2.0;

    for (int i = 0; i < ncount; ++i)
        new_row(degree_sense, degree_rhs);

    vector<double> col_coeffs{1.0, 1.0};
    double lb = 0.0;
    double ub = 1.0;
    
    for (CMR::Edge &e : graph.edges) {
        double objval = e.len;
        vector<int> ends{e.end[0], e.end[1]};

        add_col(objval, ends, col_coeffs, lb, ub);
    }
    
    tour_base = TourBasis(graph, best_data);

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

CoreLP::TourBasis::TourBasis(CMR::Graph &graph,
                             CMR::Data::BestGroup &best_data) try :
    best_tour_edges(vector<double>(graph.edge_count)),
    colstat(vector<int>(graph.edge_count, CPX_AT_LOWER)),
    rowstat(vector<int>(graph.node_count, CPX_AT_LOWER))
{
    vector<int> &int_tour_edges = best_data.best_tour_edges;
    
    for(int i = 0; i < int_tour_edges.size(); ++i)
        best_tour_edges[i] = int_tour_edges[i];

    int ncount = graph.node_count;
    vector<int> &tour_nodes = best_data.best_tour_nodes;

    for (int i = 0; i < ncount; ++i) {
        int end0 = min(tour_nodes[i], tour_nodes[(i + 1) % ncount]);
        int end1 = max(tour_nodes[i], tour_nodes[(i + 1) % ncount]);

        auto edge_it = graph.edge_lookup.find(IntPair(end0, end1));

        if (edge_it == graph.edge_lookup.end()) {
            cerr << "Edge " << end0 << ", " << end1
                 << " not found in edge hash.\n";
            throw logic_error("Graph does not contain all edges in tour.");
        }
        colstat[edge_it->second] = CPX_BASIC;
    }

    if ((ncount % 2) == 0) {
        int end0 = min(tour_nodes[ncount - 2], tour_nodes[ncount - 1]);
        int end1 = max(tour_nodes[ncount - 2], tour_nodes[ncount - 1]);

        auto edge_it = graph.edge_lookup.find(IntPair(end0, end1));

        if (edge_it == graph.edge_lookup.end()) {
            cerr << "Edge " << end0 << ", " << end1
                 << " not found in edge hash.\n";
            throw logic_error("Graph does not contain all edges in tour.");
        }

        colstat[edge_it->second] = CPX_AT_UPPER;

        end0 = min(tour_nodes[0], tour_nodes[ncount - 2]);
        end1 = max(tour_nodes[0], tour_nodes[ncount - 2]);

        edge_it = graph.edge_lookup.find(IntPair(end0, end1));
        
        if (edge_it == graph.edge_lookup.end()) {
            cerr << "Edge " << end0 << ", " << end1
                 << " not found in edge hash.\n";
            throw logic_error("Graph does not contain basis edge.");
        }

        colstat[edge_it->second] = CPX_BASIC;
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("TourBasis constructor failed.");
}

CMR::LP::PivType CoreLP::primal_pivot()
{
    runtime_error err("Problem in CoreLP::primal_pivot.");
}

}
}

