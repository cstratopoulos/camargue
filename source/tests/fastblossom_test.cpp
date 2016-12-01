#include "cc_fastblossoms.hpp"
#include "tests.hpp"
#include "datagroups.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

#include <catch.hpp>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

using std::cout;
using std::setw;
using std::vector;
using std::string;
using std::pair;

#ifdef PSEP_DO_TESTS

TEST_CASE("Fast blossosm via wrapper class",
	  "[fast2m][tiny]"){
  vector<string> tests{"blossom6", "comb9"};

  for(string &fname : tests){
    SECTION(fname){
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
      PSEP::Data::GraphGroup g_dat;
      PSEP::Data::BestGroup b_dat;
      PSEP::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
	
      REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
					      g_dat, b_dat, lp_edges,
					      s_dat));

      PSEP::cc::FastBlossom fb_sep(g_dat, b_dat, s_dat);
      REQUIRE(fb_sep.cutcall() == PSEP::CutFound::Yes);
    }
  }
}

TEST_CASE("Fast blossoms via concorde", "[fast2m]"){
  vector<string> tests{"blossom6"};

  for(string &fname : tests){
    SECTION(fname){
      string
	probfile = "problems/" + fname + ".tsp",
	solfile = "test_data/tours/" + fname + ".sol",
	subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
      PSEP::Data::GraphGroup g_dat;
      PSEP::Data::BestGroup b_dat;
      PSEP::Data::SupportGroup s_dat;
      std::vector<double> lp_edges;
	
      REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
					      g_dat, b_dat, lp_edges,
					      s_dat));
      int ncount = s_dat.G_s.node_count;
      int ecount = s_dat.G_s.edge_count;
      int cutcount = 0;
      CCtsp_lpcut_in *cc_cut = nullptr;
      CCtsp_init_lpcut_in(cc_cut);

      vector<int> cut_elist = s_dat.support_elist;
            for(int &i : cut_elist) i = b_dat.perm[i];

      REQUIRE_FALSE(CCtsp_fastblossom(&cc_cut, &cutcount, ncount, ecount,
      				      cut_elist.data(),
      				      s_dat.support_ecap.data()));
      REQUIRE(cutcount == 2);
      int i = 0;
      for(CCtsp_lpcut_in *c = cc_cut; c; c = c->next){
	int cq_count = c->cliquecount;
	cout << ++i << "th cut, clique count: " << cq_count << "\n";
	cout << "sense: " << c->sense << ", rhs: " << c->rhs << "\n";
	for(int j = 0; j < cq_count; ++j){
	  cout << "\tclique " << j << ":\n";
	  CCtsp_lpclique *clq = &(c->cliques[j]);
	  CCtsp_segment *segnodes = clq->nodes;
	  cout << "\t\tsegcount: " << clq->segcount << "\n";
	  for(int k = 0; k < clq->segcount; ++k){
	    cout << "\t\tSeg " << k << ": "
		 << segnodes[k].lo << ", " << segnodes[k].hi << "\n";
	  }
	}
	cout << "\n";
      }
      
      CCtsp_free_lpcut_in(cc_cut);
      CC_IFFREE(cc_cut, CCtsp_lpcut_in);
    }
  }
}

#endif
