#include "tests.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "simpleDP.hpp"
#include "cuts.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include <catch.hpp>

using std::cout;
using std::setprecision;
using std::endl;

using std::vector;
using std::string;
using std::pair;

#ifdef CMR_DO_TESTS

/*
TEST_CASE ("Tiny simple DP translator tests",
	   "[tiny][simpleDP]") {

  typedef pair<string, int> ProbPair;
  vector<ProbPair> tests{ProbPair("fleisA9", 1), //one dp inequality
			 ProbPair("fleisB9", 2), //two dp inequalites
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
	CMR::Data::GraphGroup g_dat;
	CMR::Data::BestGroup b_dat;
	CMR::Data::SupportGroup s_dat;
	std::vector<double> lp_edges;
	
	REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
                                                 subtourfile, g_dat, b_dat,
                                                 lp_edges, s_dat));
	cout << "Best tour: \n";
	for(int i : b_dat.best_tour_nodes) cout << i << ", ";
	cout << "\n";

	CMR::CutQueue<CMR::dominoparity> dp_q(250);
	CMR::Cut<CMR::dominoparity> dominos(g_dat, b_dat, s_dat, dp_q);

	REQUIRE(dominos.cutcall() != 1);
	CHECK(dp_q.size() == expected_num_cuts);

	CMR::CutTranslate translator(g_dat);
	
	while(!dp_q.empty()){
	  vector<int> rmatind;
	  vector<double> rmatval;
	  char sense;
	  double rhs;
	  
	  const CMR::dominoparity &dp_cut = dp_q.peek_front();
	  vector<int> &bt = b_dat.best_tour_nodes;

	  cout << "Found dp with....\n";
	  cout << "Handle: ";
	  for(int i : dp_cut.degree_nodes) cout << bt[i] << ", "; cout << "\n";
	  cout << "Nonneg edges: ";
	  for(const IntPair &e : dp_cut.nonneg_edges)
	    cout << "(" << bt[e.first] << ", " << bt[e.second] << "), ";
	  cout << "\n";
	  cout << "Simple teeth: (" << dp_cut.used_teeth.size() << " total)\n";
	  for(CMR::SimpleTooth *T : dp_cut.used_teeth){
	    cout << T->root << ", (" << T->body_start << "..." << T->body_end
		 << ") -- ";
	    CMR::CandidateTeeth::print_tooth(*T, true, b_dat.best_tour_nodes);
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
*/

SCENARIO("Wrapping simple DP separation in SimpleDP class",
         "[simpleDP]"){
  vector<string> probs {
    "lin318", "d493", "att532", "u724",
    // "dsj1000", "pr1002",
    // "d2103", "pr2392"// ,
    // "pcb3038",
    // "rl5915", "pla7397",
    // "usa13509"
  };

  for(string &fname : probs){
    string
    probfile = "problems/" + fname + ".tsp",
    solfile = "test_data/tours/" + fname + ".sol",
    subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
    CMR::Data::GraphGroup g_dat;
    CMR::Data::BestGroup b_dat;
    CMR::Data::SupportGroup s_dat;
    vector<double> lp_edges;
    CMR::Data::Instance inst;
    CMR::Data::KarpPartition kpart;

    GIVEN("A subtour polytope LP solution for " + fname){
      THEN("We can get light simple DP inequalities"){
        REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
                                                 subtourfile, g_dat, b_dat,
                                                 lp_edges, s_dat, inst));
        int ncount = g_dat.m_graph.node_count;

        REQUIRE_NOTHROW(kpart = CMR::Data::KarpPartition(ncount,
                                                         inst.ptr(), 99));
        CMR::CutTranslate translator(g_dat);
        CMR::CutQueue<CMR::dominoparity> dp_q(100);

        CMR::Sep::SimpleDP sDP(g_dat, kpart, b_dat, s_dat, dp_q);

        REQUIRE(sDP.find_cuts());
        cout << "Cut queue now has size: " << dp_q.size() << "\n";

        int primal_found = 0;
        while(!dp_q.empty()){
          vector<int> rmatind;
          vector<double> rmatval;
          char sense;
          double rhs;
	  
          const CMR::dominoparity &dp_cut = dp_q.peek_front();
          vector<int> &bt = b_dat.best_tour_nodes;
          double tour_activity, lp_activity;
          REQUIRE_FALSE(translator.get_sparse_row(dp_cut, bt, rmatind,
                                                  rmatval, sense, rhs));

          translator.get_activity(tour_activity, b_dat.best_tour_edges,
                                  rmatind, rmatval);
          translator.get_activity(lp_activity, lp_edges, rmatind, rmatval);
	  
          REQUIRE(tour_activity <= rhs);
	  
          if(tour_activity == rhs && lp_activity > rhs)
            ++primal_found;
	  
          dp_q.pop_front();
        }
        cout << "\t" << primal_found << " primal violated cuts found.\n";
      }
    }
  }
}

#endif //CMR_DO_TESTS
