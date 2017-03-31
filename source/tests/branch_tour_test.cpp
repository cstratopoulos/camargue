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

SCENARIO ("Analyzing branch constraints with BranchTourFind",
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
        int ncount = core_graph.node_count();
        const auto &adj_list = core_graph.get_adj();

        ABC::BranchTourFind bt_find(inst, b_dat, core_graph);

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
