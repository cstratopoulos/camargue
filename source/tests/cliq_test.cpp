#include "tests.hpp"

#ifdef CMR_DO_TESTS

#include "util.hpp"
#include "datagroups.hpp"
#include "cc_lpcuts.hpp"
#include "cliq.hpp"


#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

#include <cstdlib>

#include <catch.hpp>

using std::min;
using std::max;
using std::array;
using std::vector;

using std::string;
using std::cout;



static void make_tour_perm(vector<int> &tour, vector<int> &perm)
{
    for (int i = 0; i < tour.size(); ++i)
        tour[i] = i;

    std::random_shuffle(tour.begin(), tour.end());

    for (int i = 0; i < tour.size(); ++i)
        perm[tour[i]] = i;
}

SCENARIO ("Registering abstract cliques in a bank",
          "[Clique][CliqueBank]") {
    vector<int> sizes{50, 250, 500, 1000};
    for (int sz : sizes) {
        GIVEN ("A " + std::to_string(sz) +
               " node tour and a CliqueBank and some nodes") {
            vector<int> tour(sz), perm(sz);
            make_tour_perm(tour, perm);

            CMR::Sep::CliqueBank cbank(tour, perm);

            int e1 = rand() % sz, e2 = rand() % sz;
            vector<int> nodes1;
            for (int i = min(e1, e2); i <= max(e1, e2); ++i)
                nodes1.push_back(i);

            int f1 = rand() % sz, f2 = rand() % sz;
            vector<int> nodes2;
            for (int i = min(f1, f2); i <= max(f1, f2); ++i)
                nodes2.push_back(i);

            WHEN ("We add distinct nodes") {
                auto ptr1 = cbank.add_clique(nodes1);
                auto ptr2 = cbank.add_clique(nodes2);

                THEN ("Size and refcounts increase") {
                    REQUIRE(cbank.size() == 2);
                    REQUIRE(ptr1.use_count() == 2);
                    REQUIRE(ptr2.use_count() == 2);

                    AND_WHEN ("We add a duplicate") {
                        auto ptr1_copy = cbank.add_clique(nodes1);

                        THEN ("Size is unchanged but ref goes up") {
                            REQUIRE(cbank.size() == 2);
                            REQUIRE(ptr1_copy.use_count() == 3);
                            
                            AND_WHEN ("We remove a copied clique") {
                                cbank.del_clique(ptr1_copy);
                                
                                THEN ("Size is unchanged but ref goes down") {
                                    REQUIRE(cbank.size() == 2);
                                    REQUIRE(ptr1.use_count() == 2);
                                    
                                    AND_WHEN ("We remove a lone clique") {
                                        cbank.del_clique(ptr2);
                                        
                                        THEN ("Size goes down, clique is null")
                                        {
                                            REQUIRE(cbank.size() == 1);
                                            REQUIRE_FALSE(ptr2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

SCENARIO ("Testing equality and hash values of cliques",
          "[.Clique][.hash][tiny]") {
    
    GIVEN ("A tour and some nodes") {
        vector<int> tour(10);
        vector<int> perm(10);

        make_tour_perm(tour, perm);

        vector<int> nodes{1,2,7};

        WHEN ("We construct two of the same clique") {
            CMR::Sep::Clique clq1(nodes, perm), clq2(nodes, perm);

            THEN ("They are equal and hash equal") {
                REQUIRE(clq1 == clq2);
                std::size_t h1 = std::hash<CMR::Sep::Clique>{}(clq1);
                std::size_t h2 = std::hash<CMR::Sep::Clique>{}(clq2);
                REQUIRE(h1 == h2);
            }

            AND_WHEN ("We construct another clique from permuted nodes") {
                vector<int> nodes3{2,7,1};
                CMR::Sep::Clique clq3(nodes3, perm);
                THEN ("It is equal to the others") {
                    REQUIRE(clq1 == clq3);
                }
            }
        }
    }
}

SCENARIO ("Generating cliques from Concorde cliques",
          "[Clique][Sep]") {
    vector<string> probs{"pr76", "lin318", "d493", "pr1002"};
    for (string &fname : probs) {
        GIVEN("TSP instance " + fname) {
            WHEN("We find primal blossoms") {
                THEN("We can represent them with CMR cliques") {
                    string
                    probfile = "problems/" + fname + ".tsp",
                    solfile = "test_data/tours/" + fname + ".sol",
                    blossomfile = "test_data/blossom_lp/" + fname + ".2m.x",
                    subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

                    CMR::Data::GraphGroup g_dat;
                    CMR::Data::BestGroup b_dat;
                    CMR::Data::SupportGroup s_dat;
                    std::vector<double> lp_edges;
                    CMR::Sep::LPcutList cutq;
                    REQUIRE_NOTHROW(CMR::Data::make_cut_test(probfile, solfile,
                                                             subtourfile,
                                                             g_dat, b_dat,
                                                             lp_edges, s_dat));

	
                    CMR::TourGraph TG(b_dat.best_tour_edges,
                                      g_dat.core_graph.get_edges(),
                                      b_dat.perm);
                    for (int &i : s_dat.support_elist) i = b_dat.perm[i];
	
                    CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist,
                                                  s_dat.support_ecap, TG, cutq);
                    REQUIRE(fb_sep.find_cuts());

                    vector<int> &perm = b_dat.perm;
                    vector<int> &tour = b_dat.best_tour_nodes;

                    for (auto cur = cutq.begin(); cur; cur = cur->next) {
                        CCtsp_lpcut_in &cc_cut = *cur;
                        for (int i = 0; i < cc_cut.cliquecount; ++i) {
                            CCtsp_lpclique &cc_clq = cc_cut.cliques[i];

                            CMR::Sep::Clique clq(cc_clq, tour, perm, tour);

                            vector<int> cmr_nodes = clq.node_list(tour);
                            int *ar;
                            int count;

                            REQUIRE_FALSE(CCtsp_clique_to_array(&cc_clq, &ar,
                                                                &count));

                            vector<int> cc_nodes;

                            for (int i = 0; i < count; ++i)
                                cc_nodes.push_back(tour[ar[i]]);

                            std::sort(cmr_nodes.begin(), cmr_nodes.end());
                            std::sort(cc_nodes.begin(), cc_nodes.end());

                            REQUIRE(cmr_nodes == cc_nodes);
                  
                            free(ar);
                        }
                    }
                }
            }
        }
    }
}                

SCENARIO ("Abstract testing of tiny printed cliques",
          "[.Clique][.Sep][tiny]") {
    vector<int> test_sizes{10, 12, 14, 16};

    for (int sz : test_sizes) {
        GIVEN ("A tour on " + std::to_string(sz) + " nodes.") {
            vector<int> tour(sz);
            vector<int> perm(sz);

            for (int i = 0; i < sz; ++i)
                tour[i] = i;

            std::random_shuffle(tour.begin(), tour.end());

            for (int i = 0; i < sz; ++i)
                perm[tour[i]] = i;

            cout << "Tour:\n";
            for (int i : tour) cout << i << ", ";
            cout << "\n";

            cout << "Perm:\n";
            for (int i : perm) cout << i << ", ";
            cout << "\n";
            WHEN ("We grab contiguous sections") {
                vector<int> cts1(tour.begin(), tour.begin() + (sz / 4));
                vector<int> cts2(tour.begin() + (sz / 3),
                                 tour.begin() + ((2 * sz)/3));
                vector<int> cts3(tour.end() - (sz/5), tour.end());

                vector<vector<int>> cts_nodes{cts1, cts2, cts3};

                THEN ("We get a single-segment clique w the expected size") {
                    for (auto &vec : cts_nodes) {
                        vector<int> vec_copy = vec;
                        cout << "nodelist: ";
                        for(int i : vec) cout << i << ", ";
                        cout << "\n";
                        CMR::Sep::Clique clq  = CMR::Sep::Clique(vec, perm);
                        const CMR::Segment seg = clq.seg_list().front();
                        cout << "Represented as: " << seg.start << ", "
                             << seg.end << "\n";
                        
                        CHECK(clq.seg_count() == 1);
                        CHECK(seg.size() == vec_copy.size());

                        vector<int> from_clq = clq.node_list(tour);

                        CHECK(from_clq == vec_copy);
                    }
                }                
            }

            WHEN ("We grab disjoint segments") {
                vector<int> one_node{tour.front()};
                vector<int> edge{tour[2], tour[(2 * sz) / 3]};
                vector<int> nodes_and_seg{tour[0]};
                nodes_and_seg.push_back(tour[2]);
                for (auto it = tour.begin() + ((3 * sz) / 4); it != tour.end();
                     ++it)
                    nodes_and_seg.push_back(*it);

                using VecInt = std::pair<std::vector<int>, int>;

                vector<VecInt> vec_expect{VecInt(one_node, 1),
                                          VecInt(edge, 2),
                                          VecInt(nodes_and_seg, 3)};

                THEN ("We get cliques with the expected size") {
                    for (auto &vp : vec_expect) {
                        vector<int> vec_copy = vp.first;
                        int expect_sz = vp.second;

                        cout << "nodelist: ";
                        for(int i : vp.first) cout << i << ", "; cout << "\n";

                        CMR::Sep::Clique clq = CMR::Sep::Clique(vp.first, perm);
                        CHECK(clq.seg_count() == expect_sz);

                        vector<int> from_clq = clq.node_list(tour);

                        std::sort(from_clq.begin(), from_clq.end());
                        std::sort(vec_copy.begin(), vec_copy.end());

                        CHECK(from_clq == vec_copy);
                    }
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
