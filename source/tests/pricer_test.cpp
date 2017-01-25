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
        // "ulysses16",
        // "dantzig42",
        // "eil51",
        // "rat99",
        // "lin318",
        // "d493",
        // "p654",
        // "u724",
        // "dsj1000",
        "pr1002",
        // "d2103",
        // "u2319",
        // "pr2392",
        // "pcb3038"
        };

    for (string &fname : probs) {
        GIVEN ("A priceless cutting_loop run on " + fname) {
            WHEN ("We instantiate a Pricer from the final data") {
                THEN("It produces the same core edge reduced costs as CPLEX") {
                    string probfile = "problems/" + fname + ".tsp";

                    OutPrefs outprefs;        
                    Solver solver(probfile, 1485361714, outprefs);
                    
                    solver.cutting_loop(false);

                    Data::GraphGroup &g_dat =
                    const_cast<Data::GraphGroup&>(solver.graph_info());

                    LP::CoreLP &core_lp =
                    const_cast<LP::CoreLP&>(solver.get_core_lp());
                
                    Price::Pricer
                    pricer(core_lp, solver.inst_info(), g_dat);


                    vector<Price::PrEdge> pr_edges;

                    for (const Graph::Edge &e :
                         g_dat.core_graph.get_edges())
                        pr_edges.emplace_back(e.end[0], e.end[1]);

                    REQUIRE_NOTHROW(pricer.price_edges(pr_edges, true));
                    vector<double> cpx_rc =
                    core_lp.redcosts(0, core_lp.num_cols() - 1);

                    vector<double> node_pi =
                    core_lp.pi(0, solver.inst_info().node_count() - 1);
                    vector<double> cut_pi =
                    core_lp.pi(solver.inst_info().node_count(),
                               core_lp.num_rows() - 1);

                    const vector<int> &def_tour = solver.best_info().
                    best_tour_nodes;

                    for (int i = 0; i < pr_edges.size(); ++i) {
                        INFO("The edge " << i << ", ends "
                             << pr_edges[i].end[0] << ", "
                             << pr_edges[i].end[1] << ", len "
                             << g_dat.core_graph.get_edge(i).len);
                        CHECK(pr_edges[i].redcost ==
                              Approx(cpx_rc[i]).epsilon(Epsilon::Zero));
                        if (abs(pr_edges[i].redcost - cpx_rc[i]) >=
                            Epsilon::Zero) {
                            
                            int dom_t_ct = 0;
                            int dom_h_ct = 0;
                            int subct = 0;
                            int combct = 0;
                            const vector<int> &perm =
                            core_lp.external_cuts().get_cbank().ref_perm();
                            using CutType = Sep::HyperGraph::Type;

                            int end0 = pr_edges[i].end[0];
                            int end1 = pr_edges[i].end[1];
                            
                            int e0 = perm[end0];
                            int e1 = perm[end1];
                            
                            for (const Sep::HyperGraph &H :
                                 core_lp.external_cuts().get_cuts()) {
                                switch(H.cut_type()) {
                                case CutType::Subtour:
                                    if (H.get_cliques()[0]->contains(e0) !=
                                        H.get_cliques()[0]->contains(e1))
                                        ++subct;
                                    break;
                                case CutType::Comb:
                                    for (const Sep::Clique::Ptr &clq :
                                         H.get_cliques())
                                        if (clq->contains(e0) !=
                                            clq->contains(e1))
                                            ++combct;
                                    break;
                                case CutType::Domino:
                                    if (H.get_cliques()[0]->contains(e0) !=
                                        H.get_cliques()[0]->contains(e1))
                                        ++dom_h_ct;
                                    for (const Sep::Tooth::Ptr &T :
                                         H.get_teeth()) {
                                        if ((T->set_pair()[0].contains(e0) &&
                                             T->set_pair()[1].contains(e1)) || 
                                            (T->set_pair()[0].contains(e1) &&
                                             T->set_pair()[1].contains(e0)) ||
                                            (T->set_pair()[1].contains(e0) &&
                                             T->set_pair()[1].contains(e1)))
                                            ++dom_t_ct;
                                    }
                                    break;
                                }
                            }

                            cout << "\tcrosses "
                                 << subct << " subtours"
                                 << ", " << combct << " combs, "
                                 << dom_t_ct << " teeth, " << dom_h_ct
                                 << " handles.\n"
                                 << "Edge node pis: "
                                 << node_pi[end0] << ", " << node_pi[end1]
                                 << "\n";
                        }
                    }
                }
            }
        }
    }
}

/*
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
*/

#endif //CMR_DO_TESTS
