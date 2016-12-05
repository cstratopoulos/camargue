#include "tests.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "DPgraph.hpp"
#include "err_util.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <catch.hpp>


using std::cout;
using std::vector;
using std::unique_ptr;
using std::string;
using std::pair;

#ifdef PSEP_DO_TESTS

SCENARIO("Karp partition cutgraph tests",
	 "[karp][cutgraph]"){
  vector<string> probs{
    // "dsj1000", "pr1002",
    // "rl1304", "d2103"
    //"pcb3038"//, "rl5915",
    "pla7397"//,
    //"usa13509"
    //"lin318", "d493", "att532",
    //"pr1002"
  };

  cout << "Size of tooth: " << sizeof(PSEP::SimpleTooth) << "\n"
       << "Size of tooth ptr: " << sizeof(PSEP::SimpleTooth::Ptr) << "\n";
  for(string &fname : probs){
    string
      probfile = "problems/" + fname + ".tsp",
      solfile = "test_data/tours/" + fname + ".sol",
      subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
    PSEP::Data::GraphGroup g_dat;
    PSEP::Data::BestGroup b_dat;
    PSEP::Data::SupportGroup s_dat;
    vector<double> lp_edges;
    unique_ptr<PSEP::Data::Instance> inst;

    GIVEN("The TSPLIB instance " + fname){
      WHEN("The problem is loaded and we have the Instance"){
	THEN("We can get a Karp partition of the data"){
	  REQUIRE_NOTHROW(PSEP::Data::make_cut_test(probfile, solfile,
						    subtourfile, g_dat, b_dat,
						    lp_edges, s_dat, inst));
	  int partcount = 0;
	  int ncount = g_dat.m_graph.node_count;
	  int partsize = 2 * sqrt(ncount);
	  //int partsize = 0.1 * ncount;
	  CCsubdiv *subdiv_list = nullptr;
	  int **part_matrix = nullptr;
	  auto sg = PSEP::make_guard([&]{
	    CC_IFFREE(subdiv_list, CCsubdiv);
	      if(part_matrix)
		for(int i = 0; i < partcount; ++i)
		  CC_IFFREE(part_matrix[i], int);
	    CC_FREE(part_matrix, int *);
	  });
	  
	  CCrandstate rstate;
	  CCutil_sprand(99, &rstate);

	  REQUIRE_FALSE(CCutil_karp_partition(ncount, inst->ptr(), partsize,
					      &partcount, &subdiv_list,
					      &part_matrix, &rstate));

	  PSEP::CandidateTeeth cands(g_dat, b_dat, s_dat);
	      
	  REQUIRE_FALSE(cands.get_light_teeth());
	  cands.unmerged_weak_elim();
	  REQUIRE_FALSE(cands.merge_and_sort());
	      
	  int orig_sz = 0;
	  for(auto &vec : cands.light_teeth) orig_sz += vec.size();
	  cout << "Got collection of " << orig_sz << " light teeth\n";
	  
	  for(int i = 0; i < partcount; ++i){
	    AND_WHEN("We consider the " + std::to_string(i) + "th part"){
	      vector<bool> part_labels(ncount, false);
	      vector<int> &perm = b_dat.perm;
	      int *partlist = part_matrix[i];
	      
	      for(int k = 0; k < subdiv_list[i].cnt; ++k){
		part_labels[perm[partlist[k]]] = true;
	      }

	      int uncleared_lists = 0;
	      int cropped_sz = 0;
	      for(int i = 0; i < cands.light_teeth.size(); ++i){
		if(part_labels[i]){
		  ++uncleared_lists;
		  cropped_sz += cands.light_teeth[i].size();
		} else
		  cands.light_teeth[i].clear();
	      }

	      REQUIRE(uncleared_lists == subdiv_list[i].cnt);
	      
	      for(auto &vec : cands.light_teeth) cropped_sz += vec.size();
	      cout << "Number of teeth after cropping: "
		   << cropped_sz << "\n\n";

	      PSEP::DPCutGraph dp_graph(cands);
	      PSEP::CutQueue<PSEP::dominoparity> dp_q(250);
	      PSEP::CutTranslate translator(g_dat);
	      
	      THEN("We can look for DP cuts in a smaller graph"){
		double ft = PSEP::zeit();
		REQUIRE(dp_graph.simple_DP_sep(dp_q) != 1);
		ft = PSEP::zeit() - ft;
		int num_found = dp_q.size();
		int num_tight = 0, num_viol = 0;
		while(!dp_q.empty()){
		  vector<int> rmatind;
		  vector<double> rmatval;
		  char sense;
		  double rhs;
	  
		  const PSEP::dominoparity &dp_cut = dp_q.peek_front();
		  vector<int> &bt = b_dat.best_tour_nodes;
		  double tour_activity, lp_activity;

		  std::set<int> nodes_used;

		  REQUIRE_FALSE(translator.get_sparse_row(dp_cut, bt, rmatind,
							  rmatval,
							  sense, rhs));

		  translator.get_activity(tour_activity, b_dat.best_tour_edges,
					  rmatind, rmatval);
		  translator.get_activity(lp_activity, lp_edges, rmatind,
					  rmatval);
	  
		  CHECK(lp_activity > rhs);
		  REQUIRE(tour_activity <= rhs);
		  num_tight += (tour_activity == rhs);
		  num_viol += (lp_activity > rhs);
		  dp_q.pop_front();
		}
		cout << "\t"
		     << fname << ": found " << num_found << " cuts, took "
		     << ft << "s. "
		     << num_tight << " tight, " << num_viol << " violated.\n";
	      }
	    }
	  }
	}
      }
    }
  }
}

/*
TEST_CASE ("simple DP cutgraph tests", "[cutgraph]") {
  

    for(string &fname : probs){
      SECTION(fname){

	
	REQUIRE_NOTHROW(PSEP::Data::make_cut_test(probfile, solfile,
						  subtourfile, g_dat, b_dat,
						  lp_edges, s_dat));
	
	PSEP::CandidateTeeth cands(g_dat.delta, g_dat.edge_marks,
				   b_dat.best_tour_nodes, b_dat.perm,
				   s_dat.G_s, s_dat.support_elist,
				   s_dat.support_ecap);
	
	REQUIRE_FALSE(cands.get_light_teeth());
	cands.weak_elim();
      
	
	
	REQUIRE(dp_graph.simple_DP_sep(dp_q) == 0);
	cout << "\n";
      }
    }
}
*/

#endif //PSEP_DO_TESTS
