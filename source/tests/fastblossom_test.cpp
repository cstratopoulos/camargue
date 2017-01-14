#include "config.hpp"
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

#ifdef CMR_DO_TESTS

SCENARIO("Primal blossom separation by fast standard heuristics",
	 "[fast2m][Sep][CutTranslate]") {
  vector<string> probs{"pr76", "d493", "pr1002"};
  for (string &fname : probs) {
    GIVEN("TSP instance " + fname) {
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	blossomfile = "test_data/blossom_lp/" + fname + ".2m.x",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

      CMR::Data::GraphGroup g_dat;
      CMR::Data::BestGroup b_dat;
      CMR::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
      CMR::Sep::LPcutList cutq;

      WHEN("The tour is good and blossoms exit") {
	THEN("Primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));

	
	  CMR::Graph::TourGraph TG(b_dat.best_tour_edges,
                            g_dat.core_graph.get_edges(),  b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist,
                                        s_dat.support_ecap, TG, cutq);
	  REQUIRE(fb_sep.find_cuts());

          vector<int> &perm = b_dat.perm;

          CMR::Sep::CutTranslate processor(g_dat);

          for (auto cur = cutq.begin(); cur; cur = cur->next) {
              vector<int> rmatind;
              vector<double> rmatval;
              char sense;
              double rhs;
              
              processor.get_sparse_row(*cur, perm, rmatind, rmatval, sense,
                                       rhs);
              double tour_act = 0, lp_act = 0;
              
              processor.get_activity(tour_act, b_dat.best_tour_edges,
                                     rmatind, rmatval);
              processor.get_activity(lp_act, lp_edges, rmatind, rmatval);
              REQUIRE(tour_act == rhs);
              REQUIRE(lp_act < rhs);
          }
	}
      }
      AND_WHEN("The tour is good but the solution is in the blossom polytope") {
      	THEN("No primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  blossomfile, g_dat, b_dat,
						  lp_edges, s_dat));

	  CMR::Graph::TourGraph TG(b_dat.best_tour_edges,
                            g_dat.core_graph.get_edges(), b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);

	  REQUIRE_FALSE(fb_sep.find_cuts());
      	}	
      }
    }
  }
}

SCENARIO("Primal heuristic fast blossom sep in tiny instances",
	 "[fast2m][tiny][Sep][CutTranslate]") {
  vector<string> probs{"blossom6", "comb9"};
  for (string &fname : probs) {
    GIVEN("TSP instance " + fname) {
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	badsolfile = "test_data/tours/" + fname + ".bad.sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

      CMR::Data::GraphGroup g_dat;
      CMR::Data::BestGroup b_dat;
      CMR::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
      CMR::Sep::LPcutList cutq;

      WHEN("The tour is good") {
	THEN("Violated primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));

	  CMR::Graph::TourGraph TG(b_dat.best_tour_edges,g_dat.core_graph.get_edges(),
			     b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);
	  REQUIRE(fb_sep.find_cuts());

          vector<int> &perm = b_dat.perm;

          CMR::Sep::CutTranslate processor(g_dat);

          for (auto cur = cutq.begin(); cur; cur = cur->next) {
              vector<int> rmatind;
              vector<double> rmatval;
              char sense;
              double rhs;
              
              processor.get_sparse_row(*cur, perm, rmatind, rmatval, sense,
                                       rhs);
              double tour_act = 0, lp_act = 0;
              
              processor.get_activity(tour_act, b_dat.best_tour_edges,
                                     rmatind, rmatval);
              processor.get_activity(lp_act, lp_edges, rmatind, rmatval);
              REQUIRE(tour_act == rhs);
              REQUIRE(lp_act < rhs);
          }
	}
      }

      AND_WHEN("The tour is bad") {
	THEN("No primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, badsolfile,
						  subtourfile,
						  g_dat, b_dat, lp_edges,
						  s_dat));

	  CMR::Graph::TourGraph TG(b_dat.best_tour_edges, g_dat.core_graph.get_edges(),
			     b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);

	  REQUIRE_FALSE(fb_sep.find_cuts());
	}
      }
    }
  }
}

#endif
