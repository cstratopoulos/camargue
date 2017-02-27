#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "exec_branch.hpp"
#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>
#include <utility>

#include <cstdlib>

#include <catch.hpp>


using std::min;
using std::max;
using std::abs;

using std::array;
using std::vector;
using std::pair;

using std::unique_ptr;

using std::string;
using std::to_string;
using std::cout;
using std::endl;

SCENARIO ("Compressing branch tours with Cliques",
          "[ABC][Executor][branch_tour][compress_tour][expand_tour]") {
    vector<string> probs{
        "burma14",
        "ulysses16",
        "gr17",
        "gr24",
        "dantzig42",
        "pr76",
        "pr152",
        "a280",
        "lin318",
        "d493",
        "gr666",
        "pr1002",
        };

    using namespace CMR;

    for (string &prob : probs) {
    GIVEN ("An Executor for " + prob) {
        OutPrefs prefs;
        Solver solver("problems/" + prob + ".tsp",
                      "test_data/tours/" + prob + ".sol",
                      9999, prefs);
        const Graph::CoreGraph &cgraph = solver.graph_info();
        const LP::ActiveTour &act_tour = solver.active_tour();
        const Data::BestGroup &best_dat = solver.best_info();
        const Data::Instance &inst = solver.inst_info();
        LP::CoreLP &core = const_cast<LP::CoreLP &>(solver.get_core_lp());

        const vector<int> &tour_edges = best_dat.best_tour_edges;
        const vector<int> &tour_nodes = best_dat.best_tour_nodes;

        int ncount = solver.inst_info().node_count();
        int ecount = cgraph.edge_count();

        ABC::Executor exec(inst, core.get_active_tour(), best_dat,
                           cgraph, core);

        int ind1 = ecount / 4;
        EndPts e1 = cgraph.get_edge(ind1);

        int ind2 = ecount / 2;
        EndPts e2 = cgraph.get_edge(ind2);

        int ind3 = ecount / 3;
        EndPts e3 = cgraph.get_edge(ind3);

        using EndsDir = ABC::Executor::EndsDir;
        vector<EndsDir> stats{
            EndsDir(e1, ABC::dir_from_int(tour_edges[ind1])),
            EndsDir(e2, ABC::dir_from_int(tour_edges[ind2])),
            EndsDir(e3, ABC::dir_from_int(tour_edges[ind3])),
            };

    THEN ("We can observe compression/tour reuse on the best tour") {
        vector<Sep::Clique::Ptr> tcliques(3, nullptr);

        vector<EndsDir> edge_stats;
        for (int i = 0; i < 3; ++i) {
            double tval;
            vector<int> btour;
            edge_stats.push_back(stats[i]);
            REQUIRE_NOTHROW(exec.branch_tour(edge_stats, tour_nodes, btour,
                                             tval));
            REQUIRE(tval == best_dat.min_tour_value);
            REQUIRE(btour == tour_nodes);
            REQUIRE_NOTHROW(tcliques[i] = exec.compress_tour(btour));
            CAPTURE(i);
            REQUIRE(tcliques[i].use_count() == (i + 2));
            vector<int> exp_tour;
            REQUIRE_NOTHROW(exp_tour = exec.expand_tour(tcliques[i]));
            REQUIRE(exp_tour == btour);
        }

      AND_THEN ("We can compress/expand contra tours") {
            for (EndsDir &ed : stats)
                ed.second = ABC::dir_from_int(1 - static_cast<int>(ed.second));
            vector<EndsDir> edge_stats{stats[0]};
            double tval;
            vector<int> btour;
            REQUIRE_NOTHROW(exec.branch_tour(edge_stats, tour_nodes, btour,
                                             tval));

            REQUIRE(btour != tour_nodes);

            Sep::Clique::Ptr tcliq;
            REQUIRE_NOTHROW(tcliq = exec.compress_tour(btour));

            vector<int> ex_btour;
            REQUIRE_NOTHROW(ex_btour = exec.expand_tour(tcliq));

            REQUIRE(ex_btour == btour);
        }
    }
    }
    }
}


SCENARIO ("Instating a Brancher and getting problems",
          "[ABC][Executor][branch_tour][branch_edge]") {
    vector<string> probs{
        "dantzig42",
        "pr76",
        "pr152",
        "a280",
        "lin318",
        "d493",
        "gr666",
        "pr1002",
        };

    using namespace CMR;

    for (string &prob : probs) {
    GIVEN ("A cutting loop run on " + prob) {
    THEN ("We can get an edge to branch on") {
        OutPrefs prefs;
        Solver solver("problems/" + prob + ".tsp",
                      9999, prefs);
        int ncount = solver.inst_info().node_count();
        LP::PivType piv = solver.cutting_loop(ncount < 100, true, true);
        if (piv != LP::PivType::Frac)
            continue;

        LP::CoreLP &core = const_cast<LP::CoreLP &>(solver
                                                    .get_core_lp());
        const Graph::CoreGraph &cgraph = solver.graph_info();
        const LP::ActiveTour &act_tour = solver.active_tour();
        const Data::BestGroup &best_dat = solver.best_info();
        const Data::Instance &inst = solver.inst_info();

        unique_ptr<ABC::Executor> exec;

        REQUIRE_NOTHROW(util::ptr_reset(exec, inst,
                                        core.get_active_tour(),
                                        best_dat, cgraph, core));

        ABC::ScoreTuple branch_tuple;

        REQUIRE_NOTHROW(branch_tuple = exec->branch_edge());

        AND_THEN ("We can compute branch tours") {

        EndPts &branch_edge = branch_tuple.ends;
        int ind = cgraph.find_edge_ind(branch_edge.end[0],
                                       branch_edge.end[1]);
        int tour_entry = best_dat.best_tour_edges[ind];
        using EndsDir = ABC::Executor::EndsDir;
        using BDir = ABC::BranchNode::Dir;

        BDir agree_dir = tour_entry == 1 ? BDir::Up : BDir::Down;

        vector<EndsDir> agree_stats{EndsDir(branch_edge, agree_dir)};
        double agree_val;
        vector<int> agree_tour;

        cout << "Calling branch_tour to agree with best" << endl;
        REQUIRE_NOTHROW(exec->branch_tour(agree_stats,
                                          best_dat.best_tour_nodes,
                                          agree_tour,
                                          agree_val));
        REQUIRE(agree_val == best_dat.min_tour_value);
        REQUIRE(agree_tour == best_dat.best_tour_nodes);

        cout << "Calling branch_tour to disagree with best" << endl;
        BDir contra_dir = tour_entry == 1 ? BDir::Down : BDir::Up;
        vector<EndsDir> contra_stats{EndsDir(branch_edge, contra_dir)};
        double contra_val;
        vector<int> contra_tour;

        REQUIRE_NOTHROW(exec->branch_tour(contra_stats,
                                          best_dat.best_tour_nodes,
                                          contra_tour,
                                          contra_val));
        REQUIRE_FALSE(contra_tour == best_dat.best_tour_nodes);
        REQUIRE(contra_val == inst.tour_length(contra_tour));
        }
    }
    }
    }
}

#endif //CMR_DO_TESTS
