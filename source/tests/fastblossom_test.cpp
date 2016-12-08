#include "tests.hpp"
#include "cc_lpcuts.hpp"
#include "datagroups.hpp"

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
	 "[fast2m]") {
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
      CMR::Cut::LPcutList cutq;

      WHEN("The tour is good and blossoms exit") {
	THEN("Primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));

	
	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			     b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);
	  REQUIRE(fb_sep.find_cuts());
	}
      }
      AND_WHEN("The tour is good but the solution is in the blossom polytope") {
      	THEN("No primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  blossomfile, g_dat, b_dat,
						  lp_edges, s_dat));

	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			     b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);

	  REQUIRE_FALSE(fb_sep.find_cuts());
      	}	
      }
    }
  }
}

SCENARIO("Primal heuristic fast blossom sep in tiny instances",
	 "[fast2m][tiny]") {
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
      CMR::Cut::LPcutList cutq;

      WHEN("The tour is good") {
	THEN("Primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));

	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			     b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);
	  REQUIRE(fb_sep.find_cuts());
	}
      }

      AND_WHEN("The tour is bad") {
	THEN("No primal blossoms are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, badsolfile,
						  subtourfile,
						  g_dat, b_dat, lp_edges,
						  s_dat));

	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			     b_dat.perm);
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);

	  REQUIRE_FALSE(fb_sep.find_cuts());
	}
      }
    }
  }
}

/*
SCENARIO("Black box testing of tiny fast blossoms",
	  "[.fast2m][.tiny]") {
  CMR::Cut::LPcutList wrap;
  CMR::Data::GraphGroup g_dat;
  CMR::Data::BestGroup b_dat;
  CMR::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  
  vector<string> probs{"blossom6", "comb9"};
  for (string &fname : probs) {
    GIVEN("TSP instance " + fname) {
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	badsolfile = "test_data/tours/" + fname + ".bad.sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

      WHEN("The tour is good") {
	REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));

	CMR::TourGraph TG(b_dat.best_tour_edges,
			   g_dat.m_graph.edges,
			   b_dat.perm);
	for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	THEN("Odd component blossoms found, all are primal") {
	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					  s_dat.G_s.node_count,
					  s_dat.G_s.edge_count,
					  &s_dat.support_elist[0],
					  &s_dat.support_ecap[0]));
	  REQUIRE(wrap.size() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.size() == 2);
	}

	AND_THEN("Grostschel-Holland blossoms found, all are primal") {
	  REQUIRE_FALSE(CCtsp_ghfastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					    s_dat.G_s.node_count,
					    s_dat.G_s.edge_count,
					    &s_dat.support_elist[0],
					    &s_dat.support_ecap[0]));
	  REQUIRE(wrap.size() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.size() == 2);
	}
      }

      AND_WHEN("The tour is bad") {
	REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, badsolfile,
						subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));

	CMR::TourGraph TG(b_dat.best_tour_edges,
			   g_dat.m_graph.edges,
			   b_dat.perm);
	for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	THEN("Odd component blossoms found, none are primal") {
	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					  s_dat.G_s.node_count,
					  s_dat.G_s.edge_count,
					  &s_dat.support_elist[0],
					  &s_dat.support_ecap[0]));
	  REQUIRE(wrap.size() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.size() == 0);
	}

	AND_THEN("Grostschel-Holland blossoms found, none are primal") {
	  REQUIRE_FALSE(CCtsp_ghfastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					    s_dat.G_s.node_count,
					    s_dat.G_s.edge_count,
					    &s_dat.support_elist[0],
					    &s_dat.support_ecap[0]));
	  REQUIRE(wrap.size() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.size() == 0);
	}
      }
    }
  }
}
*/

#endif
