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

SCENARIO("Primal comb separation by standard block comb heuristics",
	 "[blkcomb]"){
  vector<string> probs{"pr76", "lin318", "d493", "pr1002", "d2103"};
  for(string &fname : probs){
    GIVEN("TSP instance " + fname){
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	blossomfile = "test_data/blossom_lp/" + fname + ".2m.x",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

      CMR::Data::GraphGroup g_dat;
      CMR::Data::BestGroup b_dat;
      CMR::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
      CMR::Cut::LPcutIn cutq;

      WHEN("Tour is good and solution is in the subtour polytope"){
	THEN("Primal block combs are found"){
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));

	
	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			     b_dat.perm);
	  for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::BlockCombs bc_sep(g_dat, b_dat, s_dat, TG, cutq);
	  REQUIRE(bc_sep.find_cuts());
	}
      }
    }
  }
}

SCENARIO("Primal heuristic block comb sep in tiny instances",
	 "[blkcomb][tiny]"){
  vector<string> probs{"blossom6", "comb9"};
  for(string &fname : probs){
    GIVEN("TSP instance " + fname){
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	badsolfile = "test_data/tours/" + fname + ".bad.sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

      CMR::Data::GraphGroup g_dat;
      CMR::Data::BestGroup b_dat;
      CMR::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
      CMR::Cut::LPcutIn cutq;

      WHEN("The tour is good"){
	THEN("Primal block combs are found"){
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));

	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			     b_dat.perm);
	  for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);
	  REQUIRE(fb_sep.find_cuts());
	}
      }

      WHEN("The tour is bad"){
	THEN("No primal block combs are found"){
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, badsolfile,
						   subtourfile,
						   g_dat, b_dat, lp_edges,
						   s_dat));

	  CMR::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			    b_dat.perm);
	  for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);

	  REQUIRE_FALSE(fb_sep.find_cuts());
	}
      }
    }
  }
}



#endif
