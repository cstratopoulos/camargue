#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "cc_lpcuts.hpp"
#include "datagroups.hpp"
#include "process_cuts.hpp"

#include <catch.hpp>

#include <iostream>
#include <string>
#include <utility>

using std::cout;
using std::endl;
using std::string;
using std::vector;

SCENARIO("Primal blossom separation by fast standard heuristics",
	 "[fast2m][Sep][get_row]") {
    using namespace CMR;
    vector<string> probs{"pr76", "d493", "pr1002"};
    for (string &fname : probs) {
        GIVEN("TSP instance " + fname) {
            string
            probfile = "problems/" + fname + ".tsp",
            solfile = "test_data/tours/" + fname + ".sol",
            blossomfile = "test_data/blossom_lp/" + fname + ".2m.x",
            subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

            Graph::CoreGraph core_graph;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            std::vector<double> lp_edges;
            Sep::LPcutList cutq;

            WHEN("The tour is good and blossoms exit") {
                THEN("Primal blossoms are found") {
                    REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                        subtourfile, core_graph,
                                                        b_dat, lp_edges,
                                                        s_dat));


                    Graph::TourGraph TG(b_dat.best_tour_edges,
                                        core_graph.get_edges(),
                                        b_dat.perm);
                    for (int &i : s_dat.support_elist) i = b_dat.perm[i];

                    Sep::FastBlossoms fb_sep(s_dat.support_elist,
                                                  s_dat.support_ecap, TG, cutq);
                    REQUIRE(fb_sep.find_cuts());

                    vector<int> &perm = b_dat.perm;

                    for (auto cur = cutq.begin(); cur; cur = cur->next) {
                        LP::SparseRow R = Sep::get_row(*cur, perm, core_graph);
                        double tour_act =
                        Sep::get_activity(b_dat.best_tour_edges, R);
                        double lp_act = Sep::get_activity(lp_edges, R);
                        REQUIRE(tour_act == R.rhs);
                        REQUIRE(lp_act < R.rhs);
                    }
                }
            }
            AND_WHEN("The tour is good but the solution is in the blossom polytope") {
                THEN("No primal blossoms are found") {
                    REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                        blossomfile,
                                                        core_graph, b_dat,
                                                        lp_edges, s_dat));

                    Graph::TourGraph TG(b_dat.best_tour_edges,
                                             core_graph.get_edges(), b_dat.perm);
                    for (int &i : s_dat.support_elist) i = b_dat.perm[i];

                    Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);

                    REQUIRE_FALSE(fb_sep.find_cuts());
                }
            }
        }
    }
}

SCENARIO("Primal heuristic fast blossom sep in tiny instances",
	 "[fast2m][tiny][Sep][get_row]") {
    using namespace CMR;
    vector<string> probs{"blossom6", "comb9"};
    for (string &fname : probs) {
        GIVEN("TSP instance " + fname) {
            string
            probfile = "test_data/" + fname + ".tsp",
            solfile = "test_data/tours/" + fname + ".sol",
            badsolfile = "test_data/tours/" + fname + ".bad.sol",
            subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

            Graph::CoreGraph core_graph;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            std::vector<double> lp_edges;
            Sep::LPcutList cutq;

            WHEN("The tour is good") {
                THEN("Violated primal blossoms are found") {
                    REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                        subtourfile, core_graph,
                                                        b_dat, lp_edges,
                                                        s_dat));

                    Graph::TourGraph TG(b_dat.best_tour_edges,core_graph.get_edges(),
                                             b_dat.perm);
                    for (int &i : s_dat.support_elist) i = b_dat.perm[i];

                    Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);
                    REQUIRE(fb_sep.find_cuts());

                    vector<int> &perm = b_dat.perm;


                    for (auto cur = cutq.begin(); cur; cur = cur->next) {
                        LP::SparseRow R = Sep::get_row(*cur, perm, core_graph);
                        double tour_act =
                        Sep::get_activity(b_dat.best_tour_edges, R);

                        double lp_act = Sep::get_activity(lp_edges, R);
                        REQUIRE(tour_act == R.rhs);
                        REQUIRE(lp_act < R.rhs);
                    }
                }
            }

            AND_WHEN("The tour is bad") {
                THEN("No primal blossoms are found") {
                    REQUIRE_NOTHROW(Data::make_cut_test(probfile, badsolfile,
                                                        subtourfile, core_graph,
                                                        b_dat, lp_edges,
                                                        s_dat));

                    Graph::TourGraph TG(b_dat.best_tour_edges, core_graph.get_edges(),
                                             b_dat.perm);
                    for (int &i : s_dat.support_elist) i = b_dat.perm[i];

                    Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);

                    REQUIRE_FALSE(fb_sep.find_cuts());
                }
            }
        }
    }
}

#endif
