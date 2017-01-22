#include "config.hpp"
#include "tooth.hpp"
#include "datagroups.hpp"
#include "timer.hpp"

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

using CMR::IntPair;
#ifdef CMR_DO_TESTS

static int dump_segment(double cut_val, int cut_start, int cut_end,
			void *u_data)
{
    vector<CMR::Sep::ToothBody> *vec = (vector<CMR::Sep::ToothBody> *) u_data;
    vec->emplace_back(CMR::Sep::ToothBody(cut_start, cut_end, cut_val));

  return 0;
}

TEST_CASE("New tiny candidate teeth with no elim",
	  "[.SimpleTooth][tiny]") {
    using namespace CMR; 
    vector<string> tests{"fleisA9", "fleisB9", "comb9", "ulysses16"};

    for (string &fname : tests) {
        SECTION(fname) {
            string
            probfile = "problems/" + fname + ".tsp",
            solfile = "test_data/tours/" + fname + ".sol",
            subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
            Data::GraphGroup g_dat;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            std::vector<double> lp_edges;
	
            REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile, subtourfile,
                                                     g_dat, b_dat, lp_edges,
                                                     s_dat));
            int ncount = s_dat.G_s.node_count;
      
            cout << "Best tour:\n";
            for (int i : b_dat.best_tour_nodes) cout << " " << i << ", ";
            cout << "\n";
      

            Sep::CandidateTeeth cands(g_dat, b_dat, s_dat);
            REQUIRE_NOTHROW(cands.get_light_teeth());
      
            int numfound = 0;
            for (auto &v : cands.light_teeth)
                numfound += v.size();
            cout << "\t" << numfound << " teeth found.\n";
      
            cout << "\n\n";
        }
    }
}

TEST_CASE("New candidate teeth with elim",
	  "[SimpleTooth]") {
    using namespace CMR;
    vector<string> tests{
        "lin318", "d493",
        "pr1002", "rl1304",
        "d2103", "pcb3038",
        "rl5915", "pla7397",
        "usa13509"
        };

    for (string &fname : tests) {
        SECTION(fname) {
            string
            probfile = "problems/" + fname + ".tsp",
            solfile = "test_data/tours/" + fname + ".sol",
            subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
            Data::GraphGroup g_dat;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            std::vector<double> lp_edges;
	
            REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile, subtourfile,
                                                g_dat, b_dat, lp_edges,
                                                s_dat));
            int ncount = s_dat.G_s.node_count;

            Sep::CandidateTeeth cands(g_dat, b_dat, s_dat);

            cout << "Did adj zones preprocessing.\n";
      
            REQUIRE_NOTHROW(cands.get_light_teeth());

            cands.sort_by_root();

            int numfound = 0;
            for (auto &v : cands.light_teeth) {
                numfound += v.size();
                REQUIRE(std::is_sorted(v.begin(), v.end(),
                                       [](const Sep::SimpleTooth::Ptr &S,
                                          const Sep::SimpleTooth::Ptr &T)
                                       { return S->body_size() <
                                           T->body_size();}));
            }
            cout << "Got collection of " << numfound << ".\n";


            int f_count = 0;
            int m_count = 0;
            for (auto &stat : cands.stats)
                if (stat == Sep::ListStat::Full)
                    ++f_count;
                else if (stat == Sep::ListStat::Merge)
                    ++m_count;
            cout << "Did " << f_count << " full sorts, " << m_count
                 << " merge sorts ("
                 << (ncount - (f_count + m_count)) << " untouched!)\n";

            cands.profile();
            cout << "\n";
        }
    }
}

TEST_CASE("New tiny tooth constructor with brute force tests",
	  "[adj_zones][.tiny][.SimpleTooth]") {
    using namespace CMR;
    vector<string> tests{"fleisA9", "fleisB9", "comb9", "ulysses16"};
  
    for (string &fname : tests) {
        SECTION(fname) {
            string
            probfile = "problems/" + fname + ".tsp",
            solfile = "test_data/tours/" + fname + ".sol",
            subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
            Data::GraphGroup g_dat;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            std::vector<double> lp_edges;
	
            REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile, subtourfile,
                                                g_dat, b_dat, lp_edges,
                                                s_dat));
            SupportGraph &G_s = s_dat.G_s;
            int max_deg = 0;
            Sep::CandidateTeeth cands(g_dat, b_dat, s_dat);
            int ncount = g_dat.core_graph.node_count();

            if (ncount <= 20) {
                cout << "Best tour:\n";
                for (int i : b_dat.best_tour_nodes) {
                    cout << " " << i << ", "; 
                    if (G_s.nodelist[i].s_degree > max_deg)
                        max_deg = G_s.nodelist[i].s_degree;
                }
                cout << "\n";

                for (int acc = 0; acc < max_deg; ++acc) {
                    for (int i : b_dat.best_tour_nodes) {
                        if (acc < G_s.nodelist[i].s_degree)
                            cout <<  " "
                                 << G_s.nodelist[i].adj_objs[acc].other_end
                                 << "  ";
                        else
                            cout << "    ";
                    }
                    cout << "\n";
                }
                cout << "\n";

                for (vector<int> &vec : cands.adj_zones) {
                    for (int i : vec)
                        cout << setw(2) << i << "  ";
                    cout << "\n";
                }
            }

            vector<Sep::ToothBody> seg_vec;
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
            for (auto s1 = seg_vec.begin(); s1 != seg_vec.end() - 1; ++s1) {
                for (auto s2 = s1 + 1; s2 != seg_vec.end(); ++s2) {
                    Sep::ToothBody seg1 = *s1, seg2 = *s2;
	  
                    int min_start = fmin(seg1.start, seg2.start);
                    int max_end = fmin(seg1.end, seg2.end);

                    for (int root = 0; root < ncount; ++root) {
                        if (seg1.contains(root) || seg2.contains(root))
                            continue;
	    
                        bool zone_equiv = cands.root_equivalent(root, seg1,
                                                                seg2);
                        bool brute_equiv = true;
                        int actual_vx = tour[root];
                        SNode vx = G_s.nodelist[actual_vx];

                        int d = 0, end1 = -1;
                        for (d = 0; d < vx.s_degree; ++d) {
                            end1 = perm[vx.adj_objs[d].other_end];
                            if (seg1.contains(end1) != seg2.contains(end1)) {
                                brute_equiv = false;
                                break;
                            }
                        }
	    
                        CHECK(zone_equiv == brute_equiv);
                        if (zone_equiv != brute_equiv) {
                            cout << "^^Considering root " << root
                                 << ", bodies ("
                                 << seg1.start << ", "
                                 << seg1.end << ") and (" << seg2.start
                                 << ", " << seg2.end
                                 << ")^^\n";
                            cout << "Actual: "
                                 << tour[root] << " -- ( " << tour[seg1.start]
                                 << ", " << tour[seg1.end] << ") " << "( "
                                 << tour[seg2.start] << ", "
                                 << tour[seg2.end] << ")\n";
                            cout << "Examing delta d = " << d
                                 << " of root " << root << ", actual: "
                                 << actual_vx << "\n";
                            cout << "Found disagreement with end index "
                                 << end1
                                 << ", actual: "
                                 << vx.adj_objs[d].other_end << "\n";
                            cout << "seg1 contains: "
                                 << seg1.contains(end1) << ", seg 2: "
                                 << seg2.contains(end1) << "\n";
                            IntPair s1_range = cands.get_range(root, seg1,
                                                               zones);
                            IntPair s2_range = cands.get_range(root, seg2,
                                                               zones);
                            cout << "seg1 range: " << s1_range.first << ", "
                                 << s1_range.second << " seg2 range: "
                                 << s2_range.first << ", "
                                 << s2_range.second << "\n\n";
                        }
                    }
                }
            }

        }
    }
}

#endif
