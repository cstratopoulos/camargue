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
                   const string &lp_sol_fname, Graph::CoreGraph &core_graph,
                   BestGroup &best_data, vector<double> &lp_edges,
                   SupportGroup &supp_data)
{
  Instance inst;
  make_cut_test(tsp_fname, tour_nodes_fname, lp_sol_fname, core_graph,
		best_data, lp_edges, supp_data, inst);
}

void make_cut_test(const string &tsp_fname,
                   const string &tour_nodes_fname,
                   const string &lp_sol_fname, Graph::CoreGraph &core_graph,
                   BestGroup &best_data, vector<double> &lp_edges,
                   SupportGroup &supp_data,
                   Instance &inst)
{
    runtime_error err("Problem in make_cut_test");

    try {
        inst = Instance(tsp_fname, 99);
    } CMR_CATCH_PRINT_THROW("getting Instance", err);

    int ncount = inst.node_count();

    best_data.min_tour_value = 0.0;

    try { best_data.perm.resize(ncount); }
    CMR_CATCH_PRINT_THROW("doing ncount resizes", err);

    vector<int> sup_elist;
    vector<double> sup_ecap;

    try {
        util::get_tour_nodes(ncount, best_data.best_tour_nodes,
                             tour_nodes_fname);
        util::get_lp_sol(ncount, sup_elist, sup_ecap, lp_sol_fname);
    } CMR_CATCH_PRINT_THROW("getting data from file", err);

    vector<int> &best_tour_edges = best_data.best_tour_edges;

    try {
        core_graph = Graph::CoreGraph(best_data.best_tour_nodes,
                                      inst.edgelen_func());

        best_tour_edges.resize(core_graph.edge_count(), 1);
        lp_edges.resize(core_graph.edge_count(), 0);
    } CMR_CATCH_PRINT_THROW("adding best tour edges", err);

    for (int i = 0; i < ncount; ++i)
        best_data.perm[best_data.best_tour_nodes[i]] = i;

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
            } CMR_CATCH_PRINT_THROW("pushing back new lp edge", err);
        }

        lp_edges[find_ind] = sup_ecap[i];
    }

    int e0 = best_data.best_tour_nodes[0];
    int e1 = best_data.best_tour_nodes[ncount - 2];
    int bas_ind = core_graph.find_edge_ind(e0, e1);
    if (bas_ind == -1) {
        try {
            core_graph.add_edge(e0, e1, inst.edgelen(e0, e1));
            bas_ind = core_graph.find_edge_ind(e0, e1);

            best_tour_edges.push_back(0);
            lp_edges.push_back(0);
        } CMR_CATCH_PRINT_THROW("pushing back basis edge", err);
    }

    vector<int> island;

    try {
        supp_data = Data::SupportGroup(core_graph.get_edges(),
                                       lp_edges, island,
                                       core_graph.node_count());
    } CMR_CATCH_PRINT_THROW("resetting support data", err);
}

}
}
