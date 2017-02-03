#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

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

SCENARIO ("Optimizing and computing lower bounds",
          "[Pricer][Price][exact_lb][Fixed64][primal_opt]") {
    using namespace CMR;
    vector<string> probs {
        "lin318",
        "d493",
        "pr1002",
        };

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("The solver terminates") {
                THEN ("We can optimize and compare with the exact lb") {
                    int seed = 99;
                    string probfile = "problems/" + prob + ".tsp";

                    OutPrefs outprefs;
                    Solver solver(probfile, seed, outprefs);

                    solver.cutting_loop(true, true);

                    LP::CoreLP &core = const_cast<LP::CoreLP &>(solver.
                                                                get_core_lp());
                    Data::GraphGroup &g_dat =
                    const_cast<Data::GraphGroup &>(solver.graph_info());

                    core.primal_opt();
                    auto objval = core.get_objval();

                    Price::Pricer pricer(core, solver.inst_info(), g_dat);

                    auto lb = pricer.exact_lb();
                    auto d_lb = lb.to_d();
                    
                    INFO("Primal vs dual lb\n\t"
                         << objval << "\n\t" << d_lb << "\n\t Ratio\n"
                         << (d_lb / objval));

                    REQUIRE(d_lb <= objval);
                    CHECK(d_lb == Approx(objval));
                }
            }
        }
    }    
}

SCENARIO ("Computing exact lower bounds",
          "[Pricer][Price][exact_lb][Fixed64][DualGroup]") {
    using namespace CMR;
    vector<string> probs {
        "bayg29",
        "dantzig42",
        // "pr76",
        // "lin105",
        // "a280",
        // "lin318",
        // "fl417",
        // "p654"
        };

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("The solution is optimal for its edge set") {
                int seed = 99;
                string probfile = "problems/" + prob + ".tsp";

                OutPrefs outprefs;
                Solver solver(probfile, seed, outprefs);

                solver.cutting_loop(true, true);

                LP::CoreLP &core = const_cast<LP::CoreLP &>(solver.
                                                            get_core_lp());
                Data::GraphGroup &g_dat =
                const_cast<Data::GraphGroup &>(solver.graph_info());

                Price::Pricer pricer(core, solver.inst_info(), g_dat);
                LP::DualGroup<double> dg(false, core, core.external_cuts());
                for (double &d : dg.node_pi)
                    cout << "\tNode pi " << d << "\n";
                THEN ("We can compute a valid lower bound") {
                    pricer.exact_lb();
                    
                }
            }
        }
    }
}

SCENARIO ("Comparing Pricer reduced costs to CPLEX",
          "[Price][Pricer][price_edges]") {
    using namespace CMR;
    using ProbPair = std::pair<string, int>;

vector<ProbPair> probs {
    ProbPair("ulysses16", 1485441991),
    ProbPair("eil51", 0),
    ProbPair("rat99", 0/*1485393818*/),
    ProbPair("lin318", 0/*1485393598*/),
    ProbPair("d493", 0),
    ProbPair("p654", 0),
    ProbPair("u724", 0),
    ProbPair("dsj1000", 0),
    ProbPair("pr1002", 0/*99*/),
    ProbPair("d2103", 0),
    ProbPair("pr2392", 0),
    ProbPair("pcb3038", 0),
    };

    for (ProbPair &prob : probs) {
        GIVEN ("A priceless cutting_loop run on " + prob.first) {
            WHEN ("We instantiate a Pricer from the final data") {
                THEN("It produces the same core edge reduced costs as CPLEX") {
                    string fname = prob.first;
                    int seed = prob.second;
                    string probfile = "problems/" + fname + ".tsp";

                    OutPrefs outprefs;        
                    Solver solver(probfile, seed, outprefs);
                    
                    solver.cutting_loop(false, true);

                    Data::GraphGroup &g_dat =
                    const_cast<Data::GraphGroup&>(solver.graph_info());

                    LP::CoreLP &core_lp =
                    const_cast<LP::CoreLP&>(solver.get_core_lp());

                    int ncount = g_dat.core_graph.node_count();
                    int rowcount = core_lp.num_rows();
                
                    Price::Pricer pricer(core_lp, solver.inst_info(), g_dat);
                    vector<Price::PrEdge<double>> pr_edges;
                    vector<Price::PrEdge<util::Fixed64>> f64_edges;

                    for (const Graph::Edge &e : g_dat.core_graph.get_edges()) {
                        pr_edges.emplace_back(e.end[0], e.end[1]);
                        f64_edges.emplace_back(e.end[0], e.end[1]);
                    }
                        

                    std::unique_ptr<LP::DualGroup<double>> dgp;
                    std::unique_ptr<LP::DualGroup<util::Fixed64>> f64_dgp;

                    REQUIRE_NOTHROW(pricer.price_edges(pr_edges, dgp));
                    REQUIRE_NOTHROW(pricer.price_edges(f64_edges, f64_dgp));

                    cout << "Dual feasible before getting cpx_rc: "
                         << core_lp.dual_feas() << "\n";
                    
                    vector<double> cpx_rc =
                    core_lp.redcosts(0, core_lp.num_cols() - 1);

                    vector<double> node_pi = core_lp.pi(0, ncount - 1);
                    
                    vector<double> cut_pi = (rowcount > ncount) ?
                    core_lp.pi(solver.inst_info().node_count(),
                               core_lp.num_rows() - 1) : vector<double>();

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
                        CHECK(f64_edges[i].redcost.to_d() ==
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

                    for (int i = 0; i < cut_pi.size(); ++i)
                        if (cut_pi[i] != 0.0 &&
                            core_lp.external_cuts()
                            .get_cuts()[i]
                            .cut_type() == Sep::HyperGraph::Type::Domino){
                            cout << "\t\tExample had nonzero pi domino cut.\n";
                            break;
                        }
                    cout << "\n";
                }
            }
        }
    }
}


#endif //CMR_DO_TESTS
