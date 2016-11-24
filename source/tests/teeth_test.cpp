#include "tests.hpp"
#include "n_tooth.hpp"
#include "datagroups.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

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

static double avg(int sum, int trials)
{
  return (double) sum / (double) trials;
}

static bool tooth_cmp(const PSEP::nu::SimpleTooth::Ptr &T,
		      const PSEP::nu::SimpleTooth::Ptr &S)
{
  return T->body_size() < S->body_size();
}

TEST_CASE("New tiny candidate teeth with elim",
	  "[.][tooth][tiny]") {
  vector<string> tests{"fleisA9", "fleisB9", "comb9", "ulysses16"};

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
      
      cout << "Best tour:\n";
      for(int i : b_dat.best_tour_nodes) cout << " " << i << ", ";
      cout << "\n";
      

      PSEP::nu::CandidateTeeth cands(g_dat, b_dat, s_dat);
      REQUIRE_FALSE(cands.get_light_teeth());
      int numfound = 0;

      cout << "\tLEFT ADJACENT TEETH\n";
      for(vector<PSEP::nu::SimpleTooth::Ptr> &vec : cands.right_teeth){
	numfound += vec.size();
	for(const PSEP::nu::SimpleTooth::Ptr &T : vec){
	  cands.print_tooth(*T, ncount < 20);
	}
      }

      cout << "\tRIGHT ADJACENT TEETH\n";
      for(vector<PSEP::nu::SimpleTooth::Ptr> &vec : cands.left_teeth){
	numfound += vec.size();
	for(const PSEP::nu::SimpleTooth::Ptr &T : vec){
	  cands.print_tooth(*T, ncount < 20);
	}
      }

      cout << "\tDISTANT TEETH\n";
      for(vector<PSEP::nu::SimpleTooth::Ptr> &vec : cands.dist_teeth){
	numfound += vec.size();
	for(const PSEP::nu::SimpleTooth::Ptr &T : vec){
	  cands.print_tooth(*T, ncount < 20);
	}
      }
      cout << "\n\n";
    }
  }
}

TEST_CASE("New candidate teeth with elim",
	  "[tooth]") {
  vector<string> tests{"lin318", "d493", "pr1002", "rl1304", "d2103",
		       "pcb3038"};

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

      double ct = PSEP::zeit();
      PSEP::nu::CandidateTeeth cands(g_dat, b_dat, s_dat);
      ct = PSEP::zeit() - ct;
      cout << "Did adj zones preprocessing in " << ct << "s\n";
      
      double ft = PSEP::zeit();      
      REQUIRE_FALSE(cands.get_light_teeth());
      ft = PSEP::zeit() - ft;
      
      double st = PSEP::zeit();
      int numfound = 0, msort_rval = 0;
      
      for(int i = 0; i < ncount; ++i){
	msort_rval = cands.merge_and_sort(i);
	if(msort_rval) break;
	numfound += cands.light_teeth[i].size();
      }

      REQUIRE_FALSE(msort_rval);
      
      st = PSEP::zeit() - st;
      cout << "Found " << numfound << " teeth in " << ft << "s\n";
      cout << "Sorted them in " << st << "s\n";

      int m_count = 0, s_count = 0, leftover;
      for(auto &stat : cands.stats)
	if(stat == PSEP::nu::ListStat::Merge)
	  ++m_count;
	else if(stat == PSEP::nu::ListStat::Full)
	  ++s_count;
      cout << "Did " << m_count << " merges, " << s_count << " full sorts ("
	   << (ncount - m_count - s_count) << " untouched!)\n";

      
      double et = PSEP::zeit();
      cands.weak_elim();
      et = PSEP::zeit() - et;
      int after_elim = 0;
      
      for(auto &vec : cands.light_teeth)
	after_elim += vec.size();
      
      cout << "Did weak elim on the " << (m_count + s_count)
	   << " merged/sorted roots in " << et << "s.\n";
      cout << "Got " << after_elim << " light teeth.\n";
      cout << "(eliminated "
	   << (numfound - after_elim) << ")\n";
      cout << "EVERYTHING: "
	   << (ct + ft + st + et) << "s\n";
      
      cout << "\n";
    }
  }
}

TEST_CASE("New tiny tooth constructor tests",
	  "[.][tooth][tiny]") {
  vector<string> tests{"fleisA9", "fleisB9", "comb9", "ulysses16"};
  
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
      PSEP::nu::CandidateTeeth cands(g_dat, b_dat, s_dat);
      int ncount = g_dat.m_graph.node_count;

      if(ncount <= 20){
	cout << "Best tour:\n";
	for(int i : b_dat.best_tour_nodes){
	  cout << " " << i << ", "; 
	  if(G_s.nodelist[i].s_degree > max_deg)
	    max_deg = G_s.nodelist[i].s_degree;
	}
	cout << "\n";

	for(int acc = 0; acc < max_deg; ++acc){
	  for(int i : b_dat.best_tour_nodes){
	    if(acc < G_s.nodelist[i].s_degree)
	      cout <<  " " << G_s.nodelist[i].adj_objs[acc].other_end << "  ";
	    else
	      cout << "    ";
	  }
	  cout << "\n";
	}
	cout << "\n";

	for(vector<int> &vec : cands.adj_zones){
	  for(int i : vec)
	    cout << setw(2) << i << "  ";
	  cout << "\n";
	}
      }

      vector<PSEP::nu::tooth_seg> seg_vec;
      vector<int> &tour = b_dat.best_tour_nodes;
      vector<int> &perm = b_dat.perm;
      vector<int> endmark(ncount, CC_LINSUB_BOTH_END);
      
      REQUIRE_FALSE(CCcut_linsub_allcuts(G_s.node_count, G_s.edge_count,
					 &b_dat.best_tour_nodes[0],
					 &endmark[0],
					 &s_dat.support_elist[0],
					 &s_dat.support_ecap[0],
					 2.999, &seg_vec,
					 dump_segment));

      vector<vector<int>> &zones = cands.adj_zones;
      for(auto s1 = seg_vec.begin(); s1 != seg_vec.end() - 1; ++s1){
	for(auto s2 = s1 + 1; s2 != seg_vec.end(); ++s2){
	  PSEP::nu::tooth_seg seg1 = *s1, seg2 = *s2;
	  
	  int min_start = fmin(seg1.start, seg2.start);
	  int max_end = fmin(seg1.end, seg2.end);

	  for(int root = 0; root < ncount; ++root){
	    if(seg1.contains(root) || seg2.contains(root)) continue;
	    
	    bool zone_equiv = cands.root_equivalent(root, seg1, seg2);
	    bool brute_equiv = true;
	    int actual_vx = tour[root];
	    PSEP::SNode vx = G_s.nodelist[actual_vx];

	    int d = 0, end1 = -1;
	    for(d = 0; d < vx.s_degree; ++d){
	      end1 = perm[vx.adj_objs[d].other_end];
	      if(seg1.contains(end1) != seg2.contains(end1)){
		brute_equiv = false;
		break;
	      }
	    }
	    
	    CHECK(zone_equiv == brute_equiv);
	    if(zone_equiv != brute_equiv){
	      cout << "^^Considering root " << root << ", bodies ("
		   << seg1.start << ", "
		   << seg1.end << ") and (" << seg2.start << ", " << seg2.end
		   << ")^^\n";
	      cout << "Actual: " << tour[root] << " -- ( " << tour[seg1.start]
		   << ", " << tour[seg1.end] << ") " << "( "
		   << tour[seg2.start] << ", " << tour[seg2.end] << ")\n";
	      cout << "Examing delta d = " << d
		   << " of root " << root << ", actual: "
		   << actual_vx << "\n";
	      cout << "Found disagreement with end index " << end1
		   << ", actual: " << vx.adj_objs[d].other_end << "\n";
	      cout << "seg1 contains: " << seg1.contains(end1) << ", seg 2: "
		   << seg2.contains(end1) << "\n";
	      IntPair s1_range, s2_range;
	      cands.get_range(root, seg1, s1_range, zones);
	      cands.get_range(root, seg2, s2_range, zones);
	      cout << "seg1 range: " << s1_range.first << ", "
		   << s1_range.second << " seg2 range: "
		   << s2_range.first << ", " << s2_range.second << "\n\n";
	    }
	  }
	}
      }

    }
  }
}


#endif
