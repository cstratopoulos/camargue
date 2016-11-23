#include "tests.hpp"
#include "n_tooth.hpp"
#include "datagroups.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>

#include <catch.hpp>

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

using std::cout;
using std::setw;
using std::vector;
using std::string;
using std::pair;

#ifdef PSEP_DO_TESTS

static int dump_segment(double cut_val, int cut_start, int cut_end,
			void *u_data)
{
  vector<PSEP::nu::tooth_seg> *vec = (vector<PSEP::nu::tooth_seg> *) u_data;
  vec->emplace_back(PSEP::nu::tooth_seg(cut_start, cut_end, cut_val));

  return 0;
}

TEST_CASE("New tooth constructor tests", "[tooth]") {
  vector<string> tests{"fleisA9", "fleisB9", "comb9"};
  
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
      PSEP::SupportGraph &G_s = s_dat.G_s;
      int max_deg = 0;

      cout << "Best tour:\n";
      for(int i : b_dat.best_tour_nodes){
	cout << i << ", "; 
	if(G_s.nodelist[i].s_degree > max_deg)
	  max_deg = G_s.nodelist[i].s_degree;
      }
      cout << "\n";

      for(int acc = 0; acc < max_deg; ++acc){
	for(int i : b_dat.best_tour_nodes){
	  if(acc < G_s.nodelist[i].s_degree)
	    cout << G_s.nodelist[i].adj_objs[acc].other_end << "  ";
	  else
	    cout << "   ";
	}
	cout << "\n";
      }
      cout << "\n";

      PSEP::nu::CandidateTeeth cands(g_dat, b_dat, s_dat);
      int ncount = g_dat.m_graph.node_count;

      if(ncount <= 20)
	for(vector<int> &vec : cands.adj_zones){
	  for(int i : vec)
	    cout << fabs(i) << "  ";
	  cout << "\n";
	}

      vector<PSEP::nu::tooth_seg> seg_vec;
      vector<int> &tour = b_dat.best_tour_nodes;
      vector<int> &perm = b_dat.perm;
      
      REQUIRE_FALSE(CCcut_linsub_allcuts(G_s.node_count, G_s.edge_count,
					 &b_dat.best_tour_nodes[0],
					 &cands.endmark[0],
					 &s_dat.support_elist[0],
					 &s_dat.support_ecap[0],
					 2.999, &seg_vec,
					 dump_segment));

      for(auto s1 = seg_vec.begin(); s1 != seg_vec.end() - 1; ++s1){
	for(auto s2 = s1 + 1; s2 != seg_vec.end(); ++s2){
	  PSEP::nu::tooth_seg seg1 = *s1, seg2 = *s2;
	  
	  int min_start = fmin(seg1.start, seg2.start);
	  int max_end = fmin(seg1.end, seg2.end);
	  PSEP::nu::tooth_seg blanket(min_start, max_end, 0);

	  for(int root = 0; root < ncount; ++root){
	    if(blanket.contains(root)) continue;
	    
	    bool zone_equiv = cands.root_equivalent(root, seg1, seg2);
	    bool brute_equiv = true;
	    int actual_vx = tour[root];
	    PSEP::SNode vx = G_s.nodelist[actual_vx];

	    for(int d = 0; d < vx.s_degree; ++d){
	      int end = perm[vx.adj_objs[d].other_end];
	      if(seg1.contains(end) != seg2.contains(end)){
		brute_equiv = false;
		break;
	      }
	    }
	    
	    CHECK(zone_equiv == brute_equiv);
	    if(zone_equiv != brute_equiv)
	      cout << "^^Considering root " << root << ", bodies ("
		   << seg1.start << ", "
		   << seg1.end << ") and (" << seg2.start << ", " << seg2.end
		   << ")^^\n\n";
	  }
	}
      }

    }
  }
}


#endif
