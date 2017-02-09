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
using CMR::IntPair;



SCENARIO ("Consructing a GraphGroup and BestGroup",
          "[Data][BestGroup][GraphGroup]") {
    vector<string> probs{"dantzig42", "lin318", "pr1002", "pcb3038"};

    for (string& prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            string probfile = "problems/" + prob + ".tsp";
            string solfile = "test_data/tours/" + prob + ".sol";

            CMR::Data::Instance inst(probfile, 99);
            CMR::Data::GraphGroup g_dat(inst);
            CMR::Data::BestGroup b_dat;
            
            WHEN ("An LK BestGroup is constructed") {
                THEN ("Its constructor does not throw") {
                    REQUIRE_NOTHROW(b_dat = CMR::Data::BestGroup(inst, g_dat));
                }
            }

            WHEN ("A tour file BestGroup is constructed") {
                THEN ("Its constructor does not throw.") {
                    REQUIRE_NOTHROW(b_dat = CMR::Data::BestGroup(inst, g_dat,
                                                                 solfile));
                }
            }
        }
    }
}

SCENARIO ("Verifying a BestGroup with a GraphGroup",
          "[Data][BestGroup][GraphGroup]") {
    vector<string> probs{"pr76", "d493", "d2103", "rl5915"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            string probfile = "problems/" + prob + ".tsp";
            string solfile = "test_data/tours/" + prob + ".sol";

            CMR::Data::Instance inst(probfile, 99);
            CMR::Data::GraphGroup g_dat(inst);
            int ncount = g_dat.core_graph.node_count();

            WHEN ("An LK tour is constructed") {
                THEN ("Its values can be verified by GraphGroup") {
                    CMR::Data::BestGroup b_dat(inst, g_dat);
                    double best_val = b_dat.min_tour_value;
                    double edge_val = 0;
                    double node_val = 0;

                    vector<int> &tour_nodes = b_dat.best_tour_nodes;
                    vector<int> &tour_edges = b_dat.best_tour_edges;

                    for (int i = 0; i < tour_edges.size(); ++i)
                        if (tour_edges[i])
                            edge_val += g_dat.core_graph.get_edges()[i].len;

                    bool failed_lookup = false;
                    
                    for (int i = 0; i < tour_nodes.size(); ++i) {
                        int e0 = tour_nodes[i];
                        int e1 = tour_nodes[(i + 1) % ncount];
                        int find_ind = g_dat.core_graph.find_edge_ind(e0, e1);

                        if (find_ind == -1) {
                            failed_lookup = true;
                            break;
                        }

                        node_val += g_dat.core_graph.get_edge(find_ind).len;
                    }

                    REQUIRE_FALSE(failed_lookup);
                    REQUIRE(edge_val == best_val);
                    REQUIRE(node_val == best_val);
                }
            }

            WHEN ("A file tour is constructed") {
                THEN ("Its values can be verified by GraphGroup") {
                    CMR::Data::BestGroup b_dat(inst, g_dat, solfile);
                    double best_val = b_dat.min_tour_value;
                    double edge_val = 0;
                    double node_val = 0;

                    vector<int> &tour_nodes = b_dat.best_tour_nodes;
                    vector<int> &tour_edges = b_dat.best_tour_edges;

                    for (int i = 0; i < tour_edges.size(); ++i)
                        if (tour_edges[i])
                            edge_val += g_dat.core_graph.get_edge(i).len;

                    bool failed_lookup = false;
                    
                    for (int i = 0; i < tour_nodes.size(); ++i) {
                        int e0 = tour_nodes[i];
                        int e1 = tour_nodes[(i + 1) % ncount];
                        int find_ind = g_dat.core_graph.find_edge_ind(e0, e1);

                        if (find_ind == -1) {
                            failed_lookup = true;
                            break;
                        }

                        node_val += g_dat.core_graph.get_edge(find_ind).len;
                    }

                    REQUIRE_FALSE(failed_lookup);
                    REQUIRE(edge_val == best_val);
                    REQUIRE(node_val == best_val);
                }
            }
        }
    }
}

#endif
