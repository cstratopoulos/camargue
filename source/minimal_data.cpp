#include "datagroups.hpp"
#include "err_util.hpp"
#include "io_util.hpp"
#include "util.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>

#include <cmath>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using std::vector;
using std::unique_ptr;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {
namespace Data {

void make_cut_test(const string &tsp_fname,
                   const string &tour_nodes_fname,
                   const string &lp_sol_fname, GraphGroup &graph_data,
                   BestGroup &best_data, vector<double> &lp_edges,
                   SupportGroup &supp_data)
{
  Instance inst;
  make_cut_test(tsp_fname, tour_nodes_fname, lp_sol_fname, graph_data,
		best_data, lp_edges, supp_data, inst);
}

void make_cut_test(const string &tsp_fname,
                   const string &tour_nodes_fname,
                   const string &lp_sol_fname, GraphGroup &graph_data,
                   BestGroup &best_data, vector<double> &lp_edges,
                   SupportGroup &supp_data,
                   Instance &inst)
{
    runtime_error err("Problem in make_cut_test");

    try {
        inst = Instance(tsp_fname, 99);
    } catch (const exception &e) {
        cerr << e.what() << "\n";
        throw err;
    }

    int ncount = inst.node_count();
    GraphUtils::CoreGraph &core_graph = graph_data.core_graph;

    best_data.min_tour_value = 0.0;
  
    try {
        graph_data.island.resize(ncount);
        graph_data.node_marks.resize(ncount, 0);
        best_data.perm.resize(ncount);
    } catch (const exception &e) {
        cerr << e.what() << " trying ncount resizes\n";
        throw err;
    }

    try {
        util::get_tour_nodes(ncount, best_data.best_tour_nodes,
                             tour_nodes_fname);
        util::get_lp_sol(ncount, supp_data.support_elist,
                         supp_data.support_ecap, lp_sol_fname);
    } catch (const exception &e) { cerr << e.what() << "\n"; throw err; }

    vector<int> &delta = graph_data.delta;
    vector<int> &best_tour_edges = best_data.best_tour_edges;

    try {
        core_graph = GraphUtils::CoreGraph(best_data.best_tour_nodes,
                                           inst.edgelen_func());

        delta.resize(core_graph.edge_count(), 0);
        best_tour_edges.resize(core_graph.edge_count(), 1);
        lp_edges.resize(core_graph.edge_count(), 0);
    } CMR_CATCH_PRINT_THROW("adding best tour edges", err);
  
    for (int i = 0; i < ncount; ++i)
        best_data.perm[best_data.best_tour_nodes[i]] = i;

    vector<int> &sup_elist = supp_data.support_elist;
    vector<double> &sup_ecap = supp_data.support_ecap;
    
    for (int i = 0; i < sup_ecap.size(); ++i) {
        int e0 = sup_elist[2 * i];
        int e1 = sup_elist[(2 * i) + 1];

        int find_ind = core_graph.find_edge_ind(e0, e1);
        
        if (find_ind == -1) {
            try {
                core_graph.add_edge(e0, e1, inst.edgelen(e0, e1));
                find_ind = core_graph.find_edge_ind(e0, e1);

                best_tour_edges.push_back(0);
                lp_edges.push_back(0);
                delta.push_back(0);
            } CMR_CATCH_PRINT_THROW("pushing back new lp edge", err);
        }
        
        lp_edges[find_ind] = sup_ecap[i];
    }

    try {
        supp_data.reset(core_graph.node_count(), core_graph.get_edges(),
                        lp_edges, graph_data.island);
    } CMR_CATCH_PRINT_THROW("resetting support data", err);
}

}
}
