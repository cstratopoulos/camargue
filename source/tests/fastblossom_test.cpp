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

TEST_CASE("Concorde finds violated odd component blossoms",
	  "[fast2m]"){
  PSEP::Cut::CCwrapper wrap;
  PSEP::Data::GraphGroup g_dat;
  PSEP::Data::BestGroup b_dat;
  PSEP::Data::SupportGroup s_dat;
  std::vector<double> lp_edges;
  
  SECTION("Cutcount changes appropriately"){
    vector<string> probs{"blossom6", "comb9"};
    for(string &fname : probs){
      SECTION(fname){
	string
	  probfile = "problems/" + fname + ".tsp",
	  solfile = "test_data/tours/" + fname + ".sol",
	  subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

	REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));

	  REQUIRE_FALSE(CCtsp_fastblossom(wrap.pass_ptr(), wrap.count_ptr(),
					  s_dat.G_s.node_count,
					  s_dat.G_s.edge_count,
					  &s_dat.support_elist[0],
					  &s_dat.support_ecap[0]));
	  REQUIRE(wrap.cut_count() > 0);
	  CCtsp_lpgraph L;
	  CCtsp_init_lpgraph_struct(&L);
	  auto sg = PSEP::make_guard([&L]{ CCtsp_free_lpgraph(&L); });

	  vector<int> tour_elist;
	  vector<int> &t_edges = b_dat.best_tour_edges;
	  
	  vector<PSEP::Edge> &edges = g_dat.m_graph.edges;

	  int ecount = 0;
	  for(int i = 0; i < t_edges.size(); ++i)
	    if(t_edges[i]){
	      ++ecount;
	      PSEP::Edge e = edges[i];
	      tour_elist.push_back(e.end[0]);
	      tour_elist.push_back(e.end[1]);
	    }

	  vector<double> d_tour;
	  for(int i : t_edges) d_tour.push_back(i);

	  CCtsp_build_lpgraph(&L, s_dat.G_s.node_count, ecount, &tour_elist[0],
			      (int *) NULL);
	  CCtsp_build_lpadj(&L, 0, ecount);
	  for(auto it = *wrap.pass_ptr(); it; it = it->next){
	    cout << "Slack of cut: "
		 << CCtsp_cutprice(&L, it, &d_tour[0]) << "\n";
	  }
      }
    }
  }
}

#endif
