#include "config.hpp"
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

SCENARIO("Exact primal separation of subtours",
	 "[segments]") {
  vector<string> probs{"d493", "att532", "pr1002", "rl1304", "d2103"};
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

      WHEN("The tour is good but the solution is in the subtour polytope") {
      	THEN("No segment cuts are found") {
	  REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
						  blossomfile, g_dat, b_dat,
						  lp_edges, s_dat));

	  CMR::Graph::TourGraph TG(b_dat.best_tour_edges, g_dat.core_graph.get_edges(),
			     b_dat.perm);
          
	  for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
	  CMR::Sep::SegmentCuts seg_sep(s_dat.support_elist,
                                        s_dat.support_ecap, TG, cutq);

	  REQUIRE_FALSE(seg_sep.find_cuts());
      	}	
      }
    }
  }
}

#endif
