#include "tests.hpp"
#include "karp.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "cuts.hpp"
#include "witness.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include <catch.hpp>

using std::cout;
using std::setprecision;
using std::vector;
using std::string;
using std::pair;

#ifdef CMR_DO_TESTS

SCENARIO("Finding simple DP inequalities via karp partition witnesses",
         "[karp][simpleDP][DPwitness]") {
  vector<string> probs {
    // "lin318", "d493", "att532", "u724",
    // "dsj1000", "pr1002",
    // "d2103", "pr2392"// ,
    "pcb3038",
    "rl5915", "pla7397",
    "usa13509"
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

    GIVEN("A karp partition and candidate teeth for " + fname){
      THEN("We can get simple DP inequalities in a mini cutgraph"){
        REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
                                                 subtourfile, g_dat, b_dat,
                                                 lp_edges, s_dat, inst));
        int ncount = g_dat.m_graph.node_count;

        REQUIRE_NOTHROW(kpart = CMR::Data::KarpPartition(ncount,
                                                         inst.ptr(), 99));
        CMR::CutTranslate translator(g_dat);
        
        double tt = CMR::zeit();
        CMR::CandidateTeeth cands(g_dat, b_dat, s_dat);	      
        REQUIRE_FALSE(cands.get_light_teeth());
        cands.unmerged_weak_elim();
        REQUIRE_FALSE(cands.merge_and_sort());
        tt = CMR::zeit() - tt;
	      
        int orig_sz = 0;
        for(auto &vec : cands.light_teeth) orig_sz += vec.size();
        cout << "Got collection of " << orig_sz << " light teeth in "
             << tt << "s\n";

        cout << "\tSeparating over " << kpart.num_parts() << " buckets.\n";
        int total_count = 0;
        double total_time = 0;
        
        for(int i = 0; i < kpart.num_parts(); ++i){
          double sep = CMR::zeit();
          CMR::DPwitness dpgraph(cands, kpart[i]);
          CMR::CutQueue<CMR::dominoparity> dp_q(25);
          
          REQUIRE_NOTHROW(dpgraph.simple_DP_sep(dp_q));          
          sep = CMR::zeit() - sep;
          if(dp_q.empty()){
            total_time += sep;
            continue;
          }

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

          total_count += primal_found; total_time += sep;
          if(total_count > 25){
            cout << "Breaking on part " << i << "\n"; break;
          }
        }
        cout << "\t" << total_count << " total cuts in "
             << total_time << "s\n\n";
      }
    }
  }
}

#endif //CMR_DO_TESTS
