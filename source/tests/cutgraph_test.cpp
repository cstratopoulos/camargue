#include "tests.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "DPgraph.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <catch.hpp>

using std::cout;
using std::vector;
using std::string;
using std::pair;

#ifdef PSEP_DO_TESTS

//for these tests, consider defining PSEP_DO_TESTS in DPgraph.hpp to generate
//a DOT file for the witness cutgraph which can be viewed with graphviz or
//some such
TEST_CASE ("Tiny simple DP cutgraph tests",
	   "[tiny][cutgraph]") {

  typedef pair<string, int> ProbPair;
  vector<ProbPair> tests{ProbPair("fleisA9", 1), //one dp inequality
			 ProbPair("fleisB9", 1), //one dp inequality
			 ProbPair("comb9", 2), //one blossom, one comb
			 ProbPair("ulysses16", 0)}; //no dp inequalities

  for(ProbPair &t_case : tests){
    string fname = t_case.first;
    int expected_num_cuts = t_case.second;
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
	cout << "Best tour: ";
	for(int i : b_dat.best_tour_nodes) cout << i << ", "; cout << "\n";

	  
	PSEP::CandidateTeeth cands(g_dat.delta, g_dat.edge_marks,
				   b_dat.best_tour_nodes, b_dat.perm,
				   s_dat.G_s, s_dat.support_elist,
				   s_dat.support_ecap);
	
	REQUIRE_FALSE(cands.get_light_teeth());
	cands.weak_elim();
	
	cands.print_collection();

	PSEP::DPCutGraph dp_graph(
#ifdef PSEP_DO_VIZ
				  fname,
#endif
				  cands);
	PSEP::CutQueue<PSEP::dominoparity> dp_q;

	REQUIRE(dp_graph.simple_DP_sep(dp_q) != 1);

	REQUIRE(dp_q.size() == expected_num_cuts);

	while(!dp_q.empty()){
	  PSEP::dominoparity dp = dp_q.peek_front();
	  vector<int> &bt = b_dat.best_tour_nodes;

	  cout << "Found dp with....\n";
	  cout << "Handle: ";
	  for(int i : dp.degree_nodes) cout << bt[i] << ", "; cout << "\n";
	  cout << "Nonneg edges: ";
	  for(IntPair &e : dp.nonneg_edges)
	    cout << "(" << bt[e.first] << ", " << bt[e.second] << "), ";
	  cout << "\n";
	  cout << "Simple teeth: (" << dp.used_teeth.size() << " total)\n";
	  for(PSEP::SimpleTooth *T : dp.used_teeth)
	    cands.print_tooth(*T, true);

	  dp_q.pop_front();
	}
	cout << "\n\n";
      }
    }
}

//these tests on medium sized examples should not be done with DO_VIZ. 
#ifndef PSEP_DO_VIZ
TEST_CASE ("simple DP cutgraph tests", "[cutgraph]") {
  vector<string> probs{"d2103", "pr1002", "rl1304"};

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
	
	REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));
	
	PSEP::CandidateTeeth cands(g_dat.delta, g_dat.edge_marks,
				   b_dat.best_tour_nodes, b_dat.perm,
				   s_dat.G_s, s_dat.support_elist,
				   s_dat.support_ecap);
	
	REQUIRE_FALSE(cands.get_light_teeth());
	cands.weak_elim();
      
	PSEP::DPCutGraph dp_graph(cands);
	PSEP::CutQueue<PSEP::dominoparity> dp_q;
	
	REQUIRE(dp_graph.simple_DP_sep(dp_q) == 0);
	cout << "\n";
      }
    }
}
#endif //PSEP_DO_VIZ

#endif //PSEP_TEST_TOOTH
