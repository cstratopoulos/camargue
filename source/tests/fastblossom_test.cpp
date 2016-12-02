#include "tests.hpp"
#include "cc_lpcuts.hpp"
#include "datagroups.hpp"
#include "err_util.hpp"

#include <catch.hpp>

#include <iostream>
#include <string>
#include <utility>

using std::cout;
using std::string;
using std::vector;

#ifdef PSEP_DO_TESTS

SCENARIO("Finding fast blossoms in tiny problems, filtering them by tour",
	  "[fast2m][tiny]"){
  PSEP::Cut::LPcutIn wrap;
  PSEP::Data::GraphGroup g_dat;
  PSEP::Data::BestGroup b_dat;
  PSEP::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  
  vector<string> probs{"blossom6", "comb9"};
  for(string &fname : probs){
    GIVEN("TSP instance " + fname){
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	badsolfile = "test_data/tours/" + fname + ".bad.sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

      WHEN("The tour is good"){
	REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile,
						subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));

	PSEP::TourGraph TG(b_dat.best_tour_edges,
			   g_dat.m_graph.edges,
			   b_dat.perm);
	for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	THEN("Odd component blossoms found, all are primal"){
	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					  s_dat.G_s.node_count,
					  s_dat.G_s.edge_count,
					  &s_dat.support_elist[0],
					  &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.cut_count() == 2);
	}

	AND_THEN("Grostschel-Holland blossoms found, all are primal"){
	  REQUIRE_FALSE(CCtsp_ghfastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					    s_dat.G_s.node_count,
					    s_dat.G_s.edge_count,
					    &s_dat.support_elist[0],
					    &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.cut_count() == 2);
	}
      }

      AND_WHEN("The tour is bad"){
	REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, badsolfile,
						subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));

	PSEP::TourGraph TG(b_dat.best_tour_edges,
			   g_dat.m_graph.edges,
			   b_dat.perm);
	for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	THEN("Odd component blossoms found, none are primal"){
	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					  s_dat.G_s.node_count,
					  s_dat.G_s.edge_count,
					  &s_dat.support_elist[0],
					  &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.cut_count() == 0);
	}

	AND_THEN("Grostschel-Holland blossoms found, none are primal"){
	  REQUIRE_FALSE(CCtsp_ghfastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					    s_dat.G_s.node_count,
					    s_dat.G_s.edge_count,
					    &s_dat.support_elist[0],
					    &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() == 2);
	  wrap.filter_primal(TG);
	  REQUIRE(wrap.cut_count() == 0);
	}
      }
    }
  }
}

#endif
