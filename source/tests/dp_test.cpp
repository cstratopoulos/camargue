#include "tests.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "DPgraph.hpp"

#include <vector>
#include <string>
#include <iostream>

#include <catch.hpp>

using std::cout;
using std::vector;
using std::string;

#ifdef PSEP_DO_TESTS

TEST_CASE ("Toy examples from paper",
	   "[simpleDP]") {
  vector<string> probs{"fleisA9", "fleisB9", "comb9", "ulysses16",
		       "dantzig42"};

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
	cout << "Best tour: ";
	for(int i : b_dat.best_tour_nodes) cout << i << ", "; cout << "\n";

	  
	PSEP::CandidateTeeth cands(g_dat.delta, g_dat.edge_marks,
				   b_dat.best_tour_nodes, b_dat.perm,
				   s_dat.G_s, s_dat.support_elist,
				   s_dat.support_ecap);
	
	REQUIRE_FALSE(cands.get_light_teeth());
	cands.weak_elim();
	
	for(vector<PSEP::SimpleTooth::Ptr> &vec : cands.light_teeth)
	  for(PSEP::SimpleTooth::Ptr &T : vec){
	    cout << "Root label: " << cands.print_label(*T, true) << "\n";
	  }

	PSEP::DPCutGraph dp_graph(cands.light_teeth, cands,
				  b_dat.perm, s_dat.G_s);
	dp_graph.ofname = fname;
	PSEP::CutQueue<PSEP::dominoparity> dp_q;

	REQUIRE(dp_graph.simple_DP_sep(dp_q) == 0);

	while(!dp_q.empty()){
	  PSEP::dominoparity dp = dp_q.peek_front();

	  cout << "Found dp with....\n";
	  cout << "Handle: ";
	  for(int i : dp.degree_nodes) cout << i << ", "; cout << "\n";
	  cout << "Nonneg edges: ";
	  for(IntPair &e : dp.nonneg_edges)
	    cout << "(" << e.first << ", " << e.second << "), ";
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

#endif //PSEP_TEST_TOOTH
