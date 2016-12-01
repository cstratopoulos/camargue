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

TEST_CASE("Basic member tests",
	  "[ccwrap]"){
  REQUIRE_NOTHROW(PSEP::Cut::CCwrapper wrap);
  
  PSEP::Cut::CCwrapper wrap;
  PSEP::Data::GraphGroup g_dat;
  PSEP::Data::BestGroup b_dat;
  PSEP::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  
  SECTION("Cutcount changes appropriately"){
    vector<string> probs{"blossom6", "comb9", "lin318", "pr1002"};
    for(string &fname : probs){
      SECTION(fname){
	string
	  probfile = "problems/" + fname + ".tsp",
	  solfile = "test_data/tours/" + fname + ".sol",
	  subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

	REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
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
