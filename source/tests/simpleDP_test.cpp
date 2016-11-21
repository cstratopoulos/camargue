#include "tests.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "simpleDP.hpp"
#include "cuts.hpp"

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
TEST_CASE ("Tiny simple DP translator tests",
	   "[tiny][simpleDP]") {

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

	PSEP::CutQueue<PSEP::dominoparity> dp_q(25);
	PSEP::Cut<PSEP::dominoparity> dominos(g_dat, b_dat, s_dat, dp_q);

	REQUIRE(dominos.cutcall() != 1);
	REQUIRE(dp_q.size() == expected_num_cuts);

	PSEP::CutTranslate translator(g_dat);
	
	while(!dp_q.empty()){
	  vector<int> rmatind;
	  vector<double> rmatval;
	  char sense;
	  double rhs;
	  
	  const PSEP::dominoparity &dp_cut = dp_q.peek_front();
	  vector<int> &bt = b_dat.best_tour_nodes;

	  cout << "Found dp with....\n";
	  cout << "Handle: ";
	  for(int i : dp_cut.degree_nodes) cout << bt[i] << ", "; cout << "\n";
	  cout << "Nonneg edges: ";
	  for(const IntPair &e : dp_cut.nonneg_edges)
	    cout << "(" << bt[e.first] << ", " << bt[e.second] << "), ";
	  cout << "\n";
	  cout << "Simple teeth: (" << dp_cut.used_teeth.size() << " total)\n";
	  for(PSEP::SimpleTooth *T : dp_cut.used_teeth){
	    cout << T->root << ", (" << T->body_start << "..." << T->body_end
		 << ") -- ";
	    PSEP::CandidateTeeth::print_tooth(*T, true, b_dat.best_tour_nodes);
	  }

	  REQUIRE_FALSE(translator.get_sparse_row(dp_cut,
						  b_dat.best_tour_nodes,
						  rmatind, rmatval, sense,
						  rhs));
	  double tour_activity, lp_activity;
	  translator.get_activity(tour_activity, b_dat.best_tour_edges,
				  rmatind, rmatval);
	  translator.get_activity(lp_activity, lp_edges, rmatind, rmatval);
	  
	  REQUIRE(lp_activity > rhs);
	  REQUIRE(tour_activity <= rhs);

	  dp_q.pop_front();
	}
	cout << "\n";
      }
    }
}

TEST_CASE ("simple DP cutgraph translator tests",
	   "[small][simpleDP]") {
  vector<string> probs{"st70", "gr48", "bays29", "swiss42"};

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
	bool found_nz = false;
	
	REQUIRE_FALSE(PSEP::Data::make_cut_test(probfile, solfile, subtourfile,
						g_dat, b_dat, lp_edges,
						s_dat));
	
	PSEP::CandidateTeeth cands(g_dat.delta, g_dat.edge_marks,
				   b_dat.best_tour_nodes, b_dat.perm,
				   s_dat.G_s, s_dat.support_elist,
				   s_dat.support_ecap);

	PSEP::CutQueue<PSEP::dominoparity> dp_q(25);
	PSEP::Cut<PSEP::dominoparity> dominos(g_dat, b_dat, s_dat, dp_q);
	
	REQUIRE_FALSE(dominos.cutcall());

	PSEP::CutTranslate translator(g_dat);

	
	while(!dp_q.empty()){
	  vector<int> rmatind;
	  vector<double> rmatval;
	  char sense;
	  double rhs;
	  
	  const PSEP::dominoparity &dp_cut = dp_q.peek_front();
	  vector<int> &bt = b_dat.best_tour_nodes;
	  double tour_activity, lp_activity;

	  cout << "|||||\n";
	  cout << "Considering dp with....\n";
	  cout << "Handle: ";
	  for(int i : dp_cut.degree_nodes)
	    cout << bt[i] << ", "; cout << "\n";
	  cout << "Nonneg edges: ";
	  for(const IntPair &e : dp_cut.nonneg_edges)
	    cout << "(" << bt[e.first] << ", " << bt[e.second] << "), ";
	  cout << "\n";
	  cout << "Simple teeth: (" << dp_cut.used_teeth.size()
	       << " total)\n";
	  for(PSEP::SimpleTooth *T : dp_cut.used_teeth){
	    cout << "\t" << T->root << ", (" << T->body_start
		 << "..." << T->body_end << ")\n";
	    PSEP::CandidateTeeth::print_tooth(*T, true, bt);
	    cout << "\n";
	  }
	  
	  
	  REQUIRE_FALSE(translator.get_sparse_row(dp_cut,
						  b_dat.best_tour_nodes,
						  rmatind, rmatval, sense,
						  rhs));

	  translator.get_activity(tour_activity, b_dat.best_tour_edges,
				  rmatind, rmatval);
	  translator.get_activity(lp_activity, lp_edges, rmatind, rmatval);
	  
	  REQUIRE(lp_activity > rhs);
	  REQUIRE(tour_activity <= rhs);

	  dp_q.pop_front();
	}
	cout << "\n";	
      }
    }
}

TEST_CASE ("Printless simple DP cutgraph translator tests",
	   "[simpleDP][medium]") {
  vector<string> probs{"lin318", "d493", "att532", "pr1002", "rl1304", "d2103",
		       "pcb3038"};

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

	PSEP::CutQueue<PSEP::dominoparity> dp_q(25);
	PSEP::Cut<PSEP::dominoparity> dominos(g_dat, b_dat, s_dat, dp_q);
	
	REQUIRE_FALSE(dominos.cutcall());

	PSEP::CutTranslate translator(g_dat);
	int primal_found = 0;
	while(!dp_q.empty()){
	  vector<int> rmatind;
	  vector<double> rmatval;
	  char sense;
	  double rhs;
	  
	  const PSEP::dominoparity &dp_cut = dp_q.peek_front();
	  vector<int> &bt = b_dat.best_tour_nodes;
	  double tour_activity, lp_activity;

	  REQUIRE_FALSE(translator.get_sparse_row(dp_cut, bt, rmatind, rmatval,
						  sense, rhs));

	  translator.get_activity(tour_activity, b_dat.best_tour_edges,
				  rmatind, rmatval);
	  translator.get_activity(lp_activity, lp_edges, rmatind, rmatval);
	  
	  REQUIRE(lp_activity > rhs);
	  REQUIRE(tour_activity <= rhs);
	  
	  if(tour_activity == rhs)
	    ++primal_found;
	  else
	    if(rhs - tour_activity != 1.0)
	      cout << "\tFound cut with slack " << (rhs - tour_activity)
		   << "\n";
	  
	  dp_q.pop_front();
	}
	cout << "\t" << primal_found << " primal violated cuts found\n";
	cout << "\n";	
      }
    }
}

#endif //PSEP_DO_TESTS
