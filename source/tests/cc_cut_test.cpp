#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "cc_lpcuts.hpp"
#include "datagroups.hpp"

#include <catch.hpp>

#include <iostream>
#include <string>
#include <utility>

using std::cout;
using std::string;
using std::vector;



SCENARIO("Filtering primal cuts frees and deletes cuts from list",
	 "[LPcutList][filter_primal]") {
  CMR::Graph::CoreGraph core_graph;
  CMR::Data::BestGroup b_dat;
  CMR::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  CMR::Sep::LPcutList cutq;

  GIVEN("Blossom 6 with no cuts primal") {
    WHEN("Cuts are found but none are primal") {
      THEN("Cutcount matches non-null count") {
	REQUIRE_NOTHROW(CMR::Data::make_cut_test("test_data/blossom6.tsp",
						"test_data/tours/blossom6.bad.sol",
						"test_data/subtour_lp/blossom6.sub.x",
						core_graph, b_dat, lp_edges, s_dat));
	CMR::Graph::TourGraph TG(b_dat.best_tour_edges, core_graph.get_edges(),
			   b_dat.perm);
	for (int &i : s_dat.support_elist) i = b_dat.perm[i];

	CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);

      REQUIRE_FALSE(fb_sep.find_cuts());
	int nncount = 0;
	for (auto it = cutq.begin(); it; it = it->next)
	  ++nncount;
	REQUIRE(nncount == cutq.size());
      }
    }
  }

  /*
  GIVEN("d493 with some cuts primal but not others") {
    WHEN("Cuts are found but some are not primal") {
      THEN("Cutcount matches non-null count") {
	REQUIRE_NOTHROW(CMR::Data::make_cut_test("problems/d493.tsp",
						"test_data/tours/d493.sol",
						"test_data/subtour_lp/d493.sub.x",
						core_graph, b_dat, lp_edges,
						s_dat));
	CMR::Graph::TourGraph TG(b_dat.best_tour_edges, core_graph.get_edges(),
			   b_dat.perm);
	for (int &i : s_dat.support_elist) i = b_dat.perm[i];

	CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist, s_dat.support_ecap, TG, cutq);
      REQUIRE(fb_sep.find_cuts());

	int nncount = 0;
	for (auto it = cutq.begin(); it; it = it->next)
	  ++nncount;
	REQUIRE(nncount == cutq.size());
      }
    }
  }
  */
}

// TEST_CASE("Basic member tests",
// 	  "[LPcutList]") {
//   CMR::Sep::LPcutList wrap;
//   CMR::Graph::CoreGraph core_graph;
//   CMR::Data::BestGroup b_dat;
//   CMR::Data::SupportGroup s_dat;
//   std::vector<double> lp_edges;

//   SECTION("Cutcount changes appropriately") {
//     vector<string> probs{"blossom6", "comb9", "lin318"};
//     for (string &fname : probs) {
//       SECTION(fname) {
// 	string
// 	  probfile = "problems/" + fname + ".tsp",
// 	  solfile = "test_data/tours/" + fname + ".sol",
// 	  subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

// 	REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile, subtourfile,
// 						core_graph, b_dat, lp_edges,
// 						s_dat));

// 	SECTION("Cuts found makes cut count positive") {
// 	  CMR::Data::FastBlossoms fb(s_dat.support_elist, s_dat.support_ecap, TG, wrap);
// 	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
// 					  s_dat.G_s.node_count,
// 					  s_dat.G_s.edge_count,
// 					  &s_dat.support_elist[0],
// 					  &s_dat.support_ecap[0]));
// 	  REQUIRE(wrap.size() > 0);
// 	}

// 	SECTION("No cuts means zero cut count") {
// 	  REQUIRE_FALSE(CCtsp_connect_cuts(wrap.pass_ptr(), wrap.count_ptr(),
// 					   s_dat.G_s.node_count,
// 					   s_dat.G_s.edge_count,
// 					   &s_dat.support_elist[0],
// 					   &s_dat.support_ecap[0]));
// 	  REQUIRE(wrap.size() == 0);
// 	}
//       }
//     }
//   }
// }

#endif
