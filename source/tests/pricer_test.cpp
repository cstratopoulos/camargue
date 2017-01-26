#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <cstdlib>

#include <catch.hpp>

using std::abs;
using std::min;
using std::max;
using std::array;
using std::vector;

using std::string;
using std::to_string;
using std::cout;

SCENARIO ("Comparing Pricer reduced costs to CPLEX",
          "[Price][Pricer][price_edges]") {
    using namespace CMR;
    vector<string> probs {
        "ulysses16",
        "eil51",
        "rat99", //1485393818
        "lin318", //1485393598
        "d493",
        "p654",
        "u724",
        "dsj1000",
        "pr1002",
        "d2103",
        "u2319",
        "pr2392",
        "pcb3038"
        };

    for (string &fname : probs) {
        GIVEN ("A priceless cutting_loop run on " + fname) {
            WHEN ("We instantiate a Pricer from the final data") {
                THEN("It produces the same core edge reduced costs as CPLEX") {
                    string probfile = "problems/" + fname + ".tsp";

                    OutPrefs outprefs;        
                    Solver solver(probfile, 0, outprefs);
                    
                    solver.cutting_loop(false);

                    Data::GraphGroup &g_dat =
                    const_cast<Data::GraphGroup&>(solver.graph_info());
                    int ncount = g_dat.core_graph.node_count();

                    LP::CoreLP &core_lp =
                    const_cast<LP::CoreLP&>(solver.get_core_lp());
                
                    Price::Pricer pricer(core_lp, solver.inst_info(), g_dat);
                    vector<Price::PrEdge> pr_edges;

                    for (const Graph::Edge &e : g_dat.core_graph.get_edges())
                        pr_edges.emplace_back(e.end[0], e.end[1]);

                    REQUIRE_NOTHROW(pricer.price_edges(pr_edges, true));

                    cout << "Dual feasible before getting cpx_rc: "
                         << core_lp.dual_feas() << "\n";
                    
                    vector<double> cpx_rc =
                    core_lp.redcosts(0, core_lp.num_cols() - 1);

                    vector<double> node_pi =
                    core_lp.pi(0, solver.inst_info().node_count() - 1);
                    
                    vector<double> cut_pi =
                    core_lp.pi(solver.inst_info().node_count(),
                               core_lp.num_rows() - 1);

                    vector<int> colstat = core_lp.col_stat();

                    const vector<int> &def_tour = solver.best_info().
                    best_tour_nodes;

                    for (int i = 0; i < pr_edges.size(); ++i) {
                        INFO("Edge " << i << " ("
                             << pr_edges[i].end[0] << ", "
                             << pr_edges[i].end[1] << "), len "
                             << g_dat.core_graph.get_edge(i).len
                             << ", node pis "
                             << node_pi[pr_edges[i].end[0]] << ", "
                             << node_pi[pr_edges[i].end[1]] << ", basic: "
                             << colstat[i]);
                        CHECK(pr_edges[i].redcost ==
                              Approx(cpx_rc[i]).epsilon(Epsilon::Zero));
                        if (abs(pr_edges[i].redcost - cpx_rc[i]) >=
                            Epsilon::Zero) {
                            vector<int> cmatind;
                            vector<double> cmatval;
                            core_lp.get_col(i, cmatind, cmatval);

                            for (int rownum : cmatind) {
                                if (rownum < ncount)
                                    cout << "degree eqn\n";
                                else
                                    cout << "Cut " << rownum << " "
                                         << core_lp.external_cuts()
                                    .get_cut(rownum)
                                    .cut_type() << ", pi val "
                                         << cut_pi[rownum - ncount] << "\n";
                            }
                        }
                    }
                }
            }
        }
    }
}


#endif //CMR_DO_TESTS
