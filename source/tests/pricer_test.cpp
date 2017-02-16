#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
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

SCENARIO ("Elminating edges after a run of cutting_loop",
          "[Price][exact_lb][elim]") {
    using namespace CMR;
    using f64 = util::Fixed64;
    vector<string> probs{
        "dantzig42",
        "pr76",
        "a280",
        "lin318",
        "d493",
        "pr1002",
        "pr2392",
        };

    for (string &prob : probs) {
        GIVEN ("A pure cutting loop run on " + prob + " with pricing") {
            int seed = 99;
            string probfile = "problems/" + prob + ".tsp";
            OutPrefs prefs;
            Solver solver(probfile, seed, prefs);

            auto piv = solver.cutting_loop(true, true, true);

            LP::CoreLP &core =
            const_cast<LP::CoreLP &>(solver.get_core_lp());

            THEN ("We can eliminate edges if non-optimal") {
                if (piv == LP::PivType::FathomedTour)
                    continue;

                core.primal_opt();
                auto objval = core.get_objval();
                cout << "\tPrimal opt objval " << objval << "\n";

                Graph::CoreGraph &core_graph =
                const_cast<Graph::CoreGraph &>(solver.graph_info());

                Price::Pricer pricer(core, solver.inst_info(), core_graph);

                f64 tourlen{solver.best_info().min_tour_value};
                f64 lb = pricer.exact_lb(false);
                cout << "\tExact dual lower bound " << lb << "\n";
                f64 gap{tourlen - lb};
                f64 cutoff{gap - 1};
                f64 negcutoff{0};
                util::add_mult(negcutoff, cutoff, -1);

                if (cutoff < 0) {
                    cout << "negative cutoff, do not elim\n";
                }

                cout << "\tZero cutoff " << cutoff << ", fixup cutoff "
                     << negcutoff << "\n";

                auto dg =
                util::make_unique<LP::DualGroup<f64>>(true, core,
                                                      core.
                                                      external_cuts());

                vector<Price::PrEdge<f64>> graph_edges;

                for (auto &e : core_graph.get_edges())
                    graph_edges.emplace_back(e.end[0], e.end[1]);

                pricer.price_edges(graph_edges, dg);

                int elimct = 0;
                int fixct = 0;

                for (int i = 0; i < graph_edges.size(); ++i) {
                    auto &e = graph_edges[i];
                    if (e.redcost > cutoff) {
                        ++elimct;
                        core.tighten_bound(i, 'E', 0.0);
                    } else if (e.redcost < negcutoff) {
                        ++fixct;
                        core.tighten_bound(i, 'E', 1.0);
                    }
                }
                cout << "\t" << elimct << " eliminated, " << fixct
                     << " fixed\n";
                core.primal_opt();
                cout << "New primal opt objval " << core.get_objval() << "\n\n";
            }

        }
    }
}

SCENARIO ("Computing dual bounds after run of cutting loop",
          "[LP][Price][Solver][dual][exact_lb]") {
    using namespace CMR;
    using f64 = util::Fixed64;
    vector<string> probs{
        "dantzig42",
        "pr76",
        "lin105",
        "d493",
        "pr1002",
        "pcb3038",
        };

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("Cutting loop is done running") {
                THEN ("We can compute an exact lb for the edge set") {
                    int seed = 99;
                    string probfile = "problems/" + prob + ".tsp";
                    OutPrefs prefs;
                    Solver solver(probfile, seed, prefs);

                    auto piv = solver.cutting_loop(false, false, true);

                    LP::CoreLP &core =
                    const_cast<LP::CoreLP &>(solver.get_core_lp());

                    if (piv != LP::PivType::FathomedTour) {
                        cout << "Suboptimal, optimizing\n";
                        core.primal_opt();
                    }

                    auto objval = core.get_objval();
                    cout << "\tPrimal opt objval " << objval << "\n";

                    int numrows = core.num_rows();
                    int numcols = core.num_cols();

                    auto d_pi = core.pi(0, core.num_rows() - 1);

                    vector<double> rhs;
                    core.get_rhs(rhs, 0, core.num_rows() - 1);

                    const auto &hcuts = core.external_cuts().get_cuts();
                    auto sense = core.senses(0, numrows - 1);

                    Graph::CoreGraph &core_graph =
                    const_cast<Graph::CoreGraph &>(solver.graph_info());

                    Price::Pricer pricer(core, solver.inst_info(), core_graph);
                    auto dg =
                    util::make_unique<LP::DualGroup<f64>>(true, core,
                                                          core.
                                                          external_cuts());

                    vector<Price::PrEdge<f64>> graph_edges;

                    for (auto &e : core_graph.get_edges())
                        graph_edges.emplace_back(e.end[0], e.end[1]);

                    pricer.price_edges(graph_edges, dg);

                    auto ex_rc_sum = f64{0.0};

                    for (auto &e : graph_edges)
                        if (e.redcost < 0)
                            ex_rc_sum -= e.redcost;

                    vector<f64> ex_pi(d_pi.begin(), d_pi.end());
                    auto ex_pi_sum = f64{0.0};
                    auto ex_lb = f64{0.0};
                    f64 fix_pi_sum{0.0};

                    int corcount = 0;

                    for (int i = 0; i < numrows; ++i)
                        if (sense[i] == 'G') {
                            if (ex_pi[i] < 0.0){
                                ++corcount;
                                ex_pi[i] = 0.0;
                            }
                        } else if (sense[i] == 'L') {
                            if (ex_pi[i] > 0.0) {
                                ++corcount;
                                ex_pi[i] = 0.0;
                            }
                        }
                    cout << "Corrected " << corcount << " duals\n";

                    for (int i = 0; i < rhs.size(); ++i)
                        util::add_mult(ex_pi_sum, ex_pi[i], rhs[i]);

                    ex_lb = ex_pi_sum - ex_rc_sum;
                    CHECK(ex_lb.to_d() == Approx(objval));
                }
            }
        }
    }
}


SCENARIO ("Trying to compute dual bounds for the degree LP",
          "[.LP][.Price][.price_edges][.dual]") {
    using namespace CMR;
    vector<string> probs{
        "dantzig42",
        "pr76",
        "d493",
        "dsj1000",
        "pr2392",
        };

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN("We primal opt the degree LP") {
                Data::Instance inst("problems/" + prob + ".tsp", 99);
                Graph::CoreGraph core_graph(inst);
                Data::BestGroup b_dat(inst, core_graph);
                LP::CoreLP core(core_graph, b_dat);
                core.primal_opt();
                auto objval = core.get_objval();
                THEN ("We can compute the same objval manually") {
                    auto pi = core.pi(0, core.num_rows() - 1);
                    auto rc = core.redcosts(0, core.num_cols() - 1);

                    vector<double> rhs;
                    core.get_rhs(rhs, 0, core.num_rows() - 1);

                    auto manual_val = 0.0;
                    for (auto i = 0; i < pi.size(); ++i)
                        manual_val += pi[i] * rhs[i];

                    INFO("Manual val before rc: " << manual_val);

                    auto sum_rc = 0.0;

                    for (double cost : rc)
                        if (cost < 0)
                            sum_rc -= cost;

                    INFO("Sum rc: " << sum_rc);

                    manual_val -= sum_rc;

                    REQUIRE(objval == Approx(manual_val));
                    AND_THEN ("We can get a comparable lb exactly") {
                        using f64 = util::Fixed64;
                        std::unique_ptr<LP::DualGroup<f64>> dg;

                        dg = util::make_unique<LP::DualGroup<f64>>(true, core,
                                                                   core.external_cuts());
                        auto pi_sum = f64{0.0};
                        auto i = 0;
                        for (auto &pi : dg->node_pi)
                            util::add_mult(pi_sum, pi, rhs[i++]);
                        INFO("Exact node pi before rc: " << pi_sum);

                        Price::Pricer price(core, inst, core_graph);

                        vector<Price::PrEdge<f64>> graph_edges;

                        for (auto &e : core_graph.get_edges())
                            graph_edges.emplace_back(e.end[0], e.end[1]);

                        price.price_edges(graph_edges, dg);

                        auto ex_rc_sum = f64{0.0};
                        for (auto &e : graph_edges)
                            if (e.redcost < 0)
                                ex_rc_sum -= e.redcost;
                        INFO("Exact sum rc: " << ex_rc_sum);

                        auto exact_lb = f64{0.0};
                        exact_lb = pi_sum - ex_rc_sum;
                        REQUIRE(exact_lb.to_d() <= manual_val);
                        AND_THEN ("We can get a lower lb over the full edges") {
                            int ncount = inst.node_count();

                            vector<Price::PrEdge<f64>> full_edges;

                            for (int i = 0; i < ncount; ++i)
                                for (int j = i + 1; j < ncount; ++j)
                                    full_edges.emplace_back(i, j);

                            price.price_edges(full_edges, dg);

                            auto full_rc_sum = f64{0.0};
                            for (auto &e : full_edges)
                                if (e.redcost < 0)
                                    full_rc_sum -= e.redcost;
                            INFO("Full ex sum rc: " << full_rc_sum);

                            auto full_ex_lb = f64{0.0};
                            full_ex_lb = pi_sum - full_rc_sum;
                            REQUIRE(full_ex_lb.to_d() <= manual_val);
                        }
                    }
                }
            }
        }
    }
}

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

                    solver.cutting_loop(true, true, true);

                    LP::CoreLP &core = const_cast<LP::CoreLP &>(solver.
                                                                get_core_lp());

                    Graph::CoreGraph &core_graph =
                    const_cast<Graph::CoreGraph &>(solver.graph_info());

                    core.primal_opt();
                    auto objval = core.get_objval();

                    Price::Pricer pricer(core, solver.inst_info(), core_graph);

                    auto lb = pricer.exact_lb(false);
                    auto d_lb = lb.to_d();


                    REQUIRE(d_lb <= objval);
                    CHECK(d_lb == Approx(objval));
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

                    solver.cutting_loop(false, true, true);

                    Graph::CoreGraph &core_graph =
                    const_cast<Graph::CoreGraph &>(solver.graph_info());

                    LP::CoreLP &core_lp =
                    const_cast<LP::CoreLP&>(solver.get_core_lp());

                    int ncount = core_graph.node_count();
                    int rowcount = core_lp.num_rows();

                    Price::Pricer pricer(core_lp, solver.inst_info(),
                                         core_graph);
                    vector<Price::PrEdge<double>> pr_edges;
                    vector<Price::PrEdge<util::Fixed64>> f64_edges;

                    for (const Graph::Edge &e : core_graph.get_edges()) {
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
                             << core_graph.get_edge(i).len
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
