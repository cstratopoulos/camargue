#include "config.hpp"

#ifdef CMR_DO_TESTS

#include <catch.hpp>

#include "datagroups.hpp"
#include "graph.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>

using std::min;
using std::max;
using std::cout;
using std::vector;
using std::unique_ptr;
using std::string;
using std::pair;


SCENARIO ("Constructing core LP adjacency lists",
          "[Graph][AdjList]") {
    using namespace CMR;
    vector<string> probs{"dantzig42", "lin318", "pr1002"};

    for (string& prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            string probfile = "problems/" + prob + ".tsp";
            string solfile = "test_data/tours/" + prob + ".sol";

            Data::Instance inst(probfile, 99);
            Data::GraphGroup g_dat(inst);

            Graph::AdjList full_alist;

            THEN ("We can construct an adj list from the edges") {
                REQUIRE_NOTHROW(full_alist =
                                Graph::AdjList(inst.node_count(),
                                                         g_dat.core_graph.get_edges()));

                AND_THEN ("Edge finding in the alist agrees with the edges") {
                    const vector<Graph::Edge> &edges =
                    g_dat.core_graph.get_edges();

                    for (int i = 0; i < edges.size(); ++i) {
                        auto e = edges[i];
                        auto find_ptr = full_alist.find_edge(e.end[0],
                                                             e.end[1]);
                        REQUIRE(find_ptr != nullptr);
                        REQUIRE(find_ptr->other_end == e.end[1]);
                        REQUIRE(find_ptr->edge_index == i);

                        auto rev_find = full_alist.find_edge(e.end[1],
                                                             e.end[0]);
                        REQUIRE(rev_find != nullptr);
                        REQUIRE(find_ptr->edge_index == rev_find->edge_index);
                        REQUIRE(rev_find->other_end == e.end[0]);
                    }
                }
            }

        }
    }
}

SCENARIO ("Constructing support adjacency lists",
          "[Graph][AdjList][SupportGraph]") {
    using namespace CMR;
    vector<string> probs {"pr76", "d493", "fl1577"};
    
    for (string &fname : probs) {
        string
        probfile = "problems/" + fname + ".tsp",
        solfile = "test_data/tours/" + fname + ".sol",
        subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
        Data::GraphGroup g_dat;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can construct an adj list for the support graph") {
                Data::make_cut_test(probfile, solfile, subtourfile, g_dat,
                                         b_dat, lp_edges, s_dat, inst);
                int ncount = g_dat.core_graph.node_count();

                Graph::AdjList sup_alist;
                vector<int> &sup_inds = s_dat.support_indices;
                const vector<Graph::Edge> &edges = g_dat.core_graph.get_edges();
                
                REQUIRE_NOTHROW(sup_alist =
                                Graph::AdjList(inst.node_count(),
                                                         edges, lp_edges,
                                                         sup_inds));

                AND_THEN ("We find only the edges that should be there") {
                    for (int i = 0; i < edges.size(); ++i) {
                        Graph::Edge e = edges[i];
                        auto found_ptr = sup_alist.find_edge(e.end[0],
                                                             e.end[1]);
                        if (lp_edges[i] < Epsilon::Zero)
                            REQUIRE(found_ptr == nullptr);
                        else {
                            REQUIRE(found_ptr != nullptr);
                            REQUIRE(found_ptr->edge_index == i);
                            REQUIRE(found_ptr->other_end == e.end[1]);
                            REQUIRE(found_ptr->val == lp_edges[i]);
                        }
                    }
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
