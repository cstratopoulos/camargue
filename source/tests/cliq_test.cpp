#include "tests.hpp"
#include "cliq.hpp"
#include "util.hpp"

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

#include <catch.hpp>

using std::array;
using std::vector;

using std::string;
using std::cout;

#ifdef CMR_DO_TESTS

SCENARIO ("Abstract testing of tiny printed cliques",
          "[Clique][Sep][tiny]") {
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
                        const CMR::Sep::segment seg = clq.seg_list().front();
                        cout << "Represented as: " << seg.start << ", "
                             << seg.end << "\n";
                        
                        CHECK(clq.seg_count() == 1);
                        CHECK(seg.size() == vec.size());

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
