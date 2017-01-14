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

using std::min;
using std::max;
using std::array;
using std::vector;

using std::string;
using std::to_string;
using std::cout;

SCENARIO ("Comparing Pricer reduced costs to CPLEX",
          "[Price][Pricer][price_edges]") {
    vector<string> probs {
        "dantzig42",
        "rat99",
        "lin318",
        };

    for (string &fname : probs) {
        GIVEN ("A priceless cutting_loop run on " + fname) {
            WHEN ("We instantiate a Pricer from the final data") {
                THEN("It produces the same core edge reduced costs as CPLEX") {
                    string probfile = "problems/" + fname + ".tsp";

                    CMR::OutPrefs outprefs;        
                    CMR::Solver solver(probfile, 99, outprefs);
                    
                    solver.cutting_loop(false);

                    CMR::Data::GraphGroup &g_dat =
                    const_cast<CMR::Data::GraphGroup&>(solver.graph_info());

                    CMR::LP::CoreLP &core_lp =
                    const_cast<CMR::LP::CoreLP&>(solver.get_core_lp());
                
                    CMR::Price::Pricer
                    pricer(core_lp, solver.inst_info(), g_dat);


                    vector<CMR::Price::PrEdge> pr_edges;

                    for (const CMR::Graph::Edge &e :
                         g_dat.core_graph.get_edges())
                        pr_edges.emplace_back(e.end[0], e.end[1]);

                    pricer.price_edges(pr_edges, true);
                    vector<double> cpx_rc =
                    core_lp.redcosts(0, core_lp.num_cols() - 1);

                    for (int i = 0; i < pr_edges.size(); ++i) {
                        CHECK(pr_edges[i].redcost == cpx_rc[i]);
                    }
                }
            }
        }
    }
}

SCENARIO ("Running the Solver cutting_loop on augmentable or optimal instances",
          "[Price][Pricer][add_edges]") {
    vector<string> probs{
        "pr76",
        "lin105",
        "pcb442",
        "p654",
    };

    for (string &prob : probs) {
        GIVEN (prob) {
            THEN ("We can run cutting loop and price on opt/aug tours") {
                CMR::OutPrefs prefs;
                CMR::Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   0, prefs);

                REQUIRE_NOTHROW(solver.cutting_loop(true));
            }
        }
    }
}

#endif //CMR_DO_TESTS
