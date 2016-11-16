#include "tests.hpp"
#include "datagroups.hpp"
#include "cutcall.hpp"

#include <vector>
#include <string>

#include <catch.hpp>

using std::vector;
using std::string;

#ifdef PSEP_DO_TESTS

TEST_CASE("Concorde subtour solutions are in subtour poly",
	  "[cutcontrol][in_subtour_poly]") {
  int rval = 0;
  bool result;
  PSEP::LP::Prefs prefs;

  vector<string> probs{"dantzig42", "pr1002", "lin105", "rl1304", "d2103",
		       "fleisA9", "fleisB9"};

  for(string &fname : probs){
    SECTION(fname){
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
      PSEP::Data::GraphGroup g_dat;
      PSEP::Data::BestGroup b_dat;
      PSEP::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
      rval = PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
				       g_dat, b_dat, lp_edges, s_dat);
      REQUIRE(rval == 0);

      PSEP::Data::LPGroup lp_dat(g_dat.m_graph, prefs, b_dat.perm);
      lp_dat.m_lp_edges = lp_edges;

      PSEP::CutControl controlla(g_dat, b_dat, lp_dat, s_dat, nullptr);

      rval = controlla.in_subtour_poly(result);
      REQUIRE(rval == 0);
      REQUIRE(result == true);
    }
  }
}

#endif //PSEP_DO_TESTS
