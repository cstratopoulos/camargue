#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "branch_tour.hpp"
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
using std::endl;

SCENARIO ("Computing branch tours",
          "[ABC][BranchTourFind][compute_tour]") {
    using namespace CMR;
    vector<string> probs{
        "lin318",
        "d493",
        "u724",
        "pr1002",
        "si1032",
        "d2103",
        "u2319",
        "pcb3038"
    };

    for (string &prob : probs) {
    GIVEN ("A BranchTourFind for " + prob) {
        int seed = 999;
        string probfile = "problems/" + prob + ".tsp";
        string solfile = "test_data/tours/" + prob + ".sol";
        OutPrefs prefs;
        prefs.save_tour = false;

        Solver solver(probfile, solfile, seed, prefs);
        auto piv = solver.cutting_loop(true, false, true);
        if (piv != LP::PivType::Frac)
            continue;

        const Data::Instance &inst = solver.inst_info();
        auto &b_dat = const_cast<Data::BestGroup &>(solver.best_info());
        auto &core_graph = const_cast<Graph::CoreGraph &>(solver.graph_info());
        auto &core_lp = const_cast<LP::CoreLP &>(solver.get_core_lp());

        int ncount = core_graph.node_count();
        const auto &adj_list = core_graph.get_adj();

        ABC::BranchTourFind bt_find(inst, b_dat, core_graph, core_lp);
        bt_find.verbose = true;

    WHEN ("We search for branch tours") {
        vector<ABC::EndsDir> constraints;
    THEN ("Various depths of affirmation return the best tour") {
        for (int i = 0; i < 5; ++i) {
            EndPts e = core_graph.get_edge(i);
            ABC::BranchNode::Dir dir =
            ABC::dir_from_int(b_dat.best_tour_edges[i]);
            constraints.emplace_back(e, dir);
            vector<int> btour;
            double btourval;
            bool found_tour;
            bool feas;

            bt_find.compute_tour(constraints, found_tour, feas, btour,
                                 btourval, false);
            REQUIRE(found_tour);
            REQUIRE(feas);
            REQUIRE(btourval == b_dat.min_tour_value);
            REQUIRE(btour == b_dat.best_tour_nodes);
        }
    AND_THEN("We can compute actual branch tours and compare qualities") {
        for (int i : {ncount / 4, ncount / 3, ncount / 2, 3 * ncount / 4}) {
            EndPts e = core_graph.get_edge(i);
            INFO("Contra-fixing edge " << e << ", ind " << i);
            ABC::BranchNode::Dir dir =
            ABC::dir_from_int(1 - b_dat.best_tour_edges[i]);
            constraints.emplace_back(e, dir);

            vector<int> sp_btour;
            double sp_tourval;
            bool sp_feas;
            bool sp_found;

            vector<int> real_btour;
            double real_tourval;
            bool real_feas;
            bool real_found;

            bt_find.compute_tour(constraints, sp_found, sp_feas, sp_btour,
                                 sp_tourval, false);
            bt_find.compute_tour(constraints, real_found, real_feas, real_btour,
                                 real_tourval, true);

            REQUIRE(sp_feas == real_feas);
            REQUIRE((int) sp_found <= (int) real_found);
            REQUIRE(real_tourval <= sp_tourval);
            if (!sp_found && !real_found) {
                cout << "Neither approach found a tour, breaking" << endl;
                break;
            }
            if (real_tourval < sp_tourval)
                cout << "Real improves on estimate!!" << endl;
        }

    }
    }
    }

    }
    }
}

SCENARIO ("Analyzing branch constraints",
         "[ABC][BranchTourFind][common_constraints][branch_constraints]"
          "[tour_compliant][obvious_infeas]")
{
    using namespace CMR;
    vector<string> probs{
        "pr76",
        "lin318",
        "gr431",
        "u724",
        "rl1304",
        };

    for (string &prob : probs) {
    GIVEN ("A BranchTourFind for " + prob) {
        int seed = 999;
        string probfile = "problems/" + prob + ".tsp";
        string solfile = "test_data/tours/" + prob + ".sol";
        OutPrefs prefs;
        prefs.save_tour = false;

        Solver solver(probfile, solfile, seed, prefs);
        auto piv = solver.cutting_loop(true, false, true);
        if (piv != LP::PivType::Frac)
            continue;

        const Data::Instance &inst = solver.inst_info();
        auto &b_dat = const_cast<Data::BestGroup &>(solver.best_info());
        auto &core_graph = const_cast<Graph::CoreGraph &>(solver.graph_info());
        auto &core_lp = const_cast<LP::CoreLP &>(solver.get_core_lp());

        int ncount = core_graph.node_count();
        const auto &adj_list = core_graph.get_adj();

        ABC::BranchTourFind bt_find(inst, b_dat, core_graph, core_lp);
        bt_find.verbose = true;

    WHEN ("There is a list of constraints") {
        vector<ABC::EndsDir> constraints;
        THEN ("We can detect obvious overfixing infeasibilities") {
            bool found_node = false;
            for (int i = 0; i < ncount; ++i) {
                const Graph::Node &n = adj_list.nodelist[i];
                if (n.degree() > 2) {
                    found_node = true;

                    for (const Graph::AdjObj &a : n.neighbors) {
                        EndPts e(i, a.other_end);
                        constraints.emplace_back(e, ABC::BranchNode::Dir::Up);
                    }

                    break;
                }
            }

            REQUIRE(found_node);
            REQUIRE(bt_find.obvious_infeas(constraints));
            AND_THEN ("Constraints without obvious infeasibilities look fine") {
                constraints.clear();
                REQUIRE(constraints.empty());
                REQUIRE_FALSE(bt_find.obvious_infeas(constraints)); //empty list

                auto &core_lp = const_cast<LP::CoreLP &>(solver.get_core_lp());

                vector<double> lp_vec = core_lp.lp_vec();
                for (int i = 0; i < lp_vec.size(); ++i) {
                    double val = lp_vec[i];
                    if (val < Epsilon::Zero) {
                        constraints.emplace_back(core_graph.get_edge(i),
                                                 ABC::BranchNode::Dir::Down);
                    } else if (val > 1 - Epsilon::Zero) {
                        constraints.emplace_back(core_graph.get_edge(i),
                                                 ABC::BranchNode::Dir::Up);
                    }
                }

                REQUIRE_FALSE(bt_find.obvious_infeas(constraints));
            }
        }

    THEN ("Tours are compliant with fix-up affirming branches") {
        const vector<int> &tour_nodes = b_dat.best_tour_nodes;
        const vector<int> &tour_edges = b_dat.best_tour_edges;

        for (int i = 0; i < tour_edges.size(); ++i)
            if (tour_edges[i] == 1)
                constraints.emplace_back(core_graph.get_edge(i),
                                         ABC::BranchNode::Dir::Up);

        REQUIRE(bt_find.tour_compliant(tour_nodes, constraints));
        AND_THEN("They remain compliant when fixing down non-tour edges") {
            for (int i = 0; i < tour_edges.size(); ++i)
                if (tour_edges[i] == 0)
                    constraints.emplace_back(core_graph.get_edge(i),
                                             ABC::BranchNode::Dir::Down);
            REQUIRE(bt_find.tour_compliant(tour_nodes, constraints));
            AND_THEN ("They are non-compliant when the fixings are inverted") {
                for (ABC::EndsDir &ed : constraints) {
                    ed.second = ABC::dir_from_int(1 -
                                                  static_cast<int>(ed.second));
                }
                REQUIRE_FALSE(bt_find.tour_compliant(tour_nodes, constraints));
            }
        }
    }
    }
    }
    }
}

#endif
