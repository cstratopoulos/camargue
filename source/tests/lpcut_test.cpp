#include "tests.hpp"
#include "cc_lpcuts.hpp"
#include "datagroups.hpp"

#include <catch.hpp>

#include <iostream>
#include <string>
#include <utility>

using std::cout;
using std::string;
using std::vector;

#ifdef PSEP_DO_TESTS

SCENARIO("Filtering primal cuts frees and deletes cuts from list",
	 "[LPcutIn][filter_primal]"){
  PSEP::Data::GraphGroup g_dat;
  PSEP::Data::BestGroup b_dat;
  PSEP::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  PSEP::Cut::LPcutIn cutq;

  GIVEN("Blossom 6 with no cuts primal"){
    WHEN("Cuts are found but none are primal"){
      THEN("Cutcount matches non-null count"){
	REQUIRE_NOTHROW(PSEP::Data::make_cut_test("problems/blossom6.tsp",
						"test_data/tours/blossom6.bad.sol",
						"test_data/subtour_lp/blossom6.sub.x",
						g_dat, b_dat, lp_edges, s_dat));
	PSEP::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			   b_dat.perm);
	for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	PSEP::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);

      REQUIRE_FALSE(fb_sep.find_cuts());
	int nncount = 0;
	for(auto it = cutq.begin(); it; it = it->next)
	  ++nncount;
	REQUIRE(nncount == cutq.cut_count());
      }
    }
  }

  GIVEN("d493 with some cuts primal but not others"){
    WHEN("Cuts are found but some are not primal"){
      THEN("Cutcount matches non-null count"){
	REQUIRE_NOTHROW(PSEP::Data::make_cut_test("problems/d493.tsp",
						"test_data/tours/d493.sol",
						"test_data/subtour_lp/d493.sub.x",
						g_dat, b_dat, lp_edges,
						s_dat));
	PSEP::TourGraph TG(b_dat.best_tour_edges, g_dat.m_graph.edges,
			   b_dat.perm);
	for(int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	PSEP::Cut::FastBlossoms fb_sep(g_dat, b_dat, s_dat, TG, cutq);
      REQUIRE(fb_sep.find_cuts());
      
	int nncount = 0;
	for(auto it = cutq.begin(); it; it = it->next)
	  ++nncount;
	REQUIRE(nncount == cutq.cut_count());
      }
    }
  }
}

TEST_CASE("Basic member tests",
	  "[LPcutIn]"){
  PSEP::Cut::LPcutIn wrap;
  PSEP::Data::GraphGroup g_dat;
  PSEP::Data::BestGroup b_dat;
  PSEP::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  
  SECTION("Cutcount changes appropriately"){
    vector<string> probs{"blossom6", "comb9", "lin318"};
    for(string &fname : probs){
      SECTION(fname){
	string
	  probfile = "problems/" + fname + ".tsp",
	  solfile = "test_data/tours/" + fname + ".sol",
	  subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

	REQUIRE_NOTHROW(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));

	SECTION("Cuts found makes cut count positive"){    
	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					  s_dat.G_s.node_count,
					  s_dat.G_s.edge_count,
					  &s_dat.support_elist[0],
					  &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() > 0);
	}

	SECTION("No cuts means zero cut count"){
	  REQUIRE_FALSE(CCtsp_connect_cuts(wrap.pass_ptr(), wrap.count_ptr(),
					   s_dat.G_s.node_count,
					   s_dat.G_s.edge_count,
					   &s_dat.support_elist[0],
					   &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() == 0);	  
	}
      }
    }
  }
}

#endif
