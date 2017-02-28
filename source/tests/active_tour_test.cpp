#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "active_tour.hpp"
#include "datagroups.hpp"
#include "core_lp.hpp"
#include "lp_interface.hpp"
#include "io_util.hpp"
#include "util.hpp"
#include "err_util.hpp"

#include <catch.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;


SCENARIO ("Constructing ActiveTours",
          "[LP][ActiveTour]") {
    using namespace CMR;
    vector<string> probs{
        "swiss42",
        "eil76",
        "rat99",
        "gr431",
        "d657",
        "pr1002",
        "rl1304",
        };

    for (string &prob : probs) {
    GIVEN ("An Instance/DataGroups/CoreLP for " + prob) {
        Data::Instance inst("problems/" + prob + ".tsp", 99);
        Graph::CoreGraph graph(inst);
        Data::BestGroup b_dat(inst, graph);
        LP::CoreLP core(graph, b_dat);

        int ncount = graph.node_count();
        int ecount = graph.edge_count();
        THEN ("We can construct and verify a Padberg-Hong basis") {
            LP::ActiveTour T;

            REQUIRE_NOTHROW(T = LP::ActiveTour(graph, b_dat));

            const LP::Basis &tbase = T.base();
            const vector<int> &tour_nodes = T.nodes();

            // have we checked the basis status for this edge.
            vector<bool> checked_basis(ecount, false);

            REQUIRE(tbase.colstat.size() == ecount);
            REQUIRE(tbase.rowstat.size() == ncount);

            if (ncount % 2 == 0) {
                int discard = graph.find_edge_ind(tour_nodes[ncount - 2],
                                                  tour_nodes[ncount - 1]);
                REQUIRE(tbase.colstat[discard] == LP::BStat::AtUpper);
                checked_basis[discard] = true;

                int newbase = graph.find_edge_ind(tour_nodes[0],
                                                  tour_nodes[ncount - 2]);
                REQUIRE(tbase.colstat[newbase] == LP::BStat::Basic);
                checked_basis[newbase] = true;
                cout << "\tChecked even nodecount basis edges" << endl;
            }

            for (int i = 0; i < ncount; ++i) {
                int ind = graph.find_edge_ind(tour_nodes[i],
                                              tour_nodes[(i + 1) % ncount]);
                if (!checked_basis[ind]) {
                    REQUIRE(tbase.colstat[ind] == LP::BStat::Basic);
                    checked_basis[ind] = true;
                }
            }
            cout << "\tChecked all tour edges" << endl;

            for (int i = 0; i < ecount; ++i)
                if (!checked_basis[i])
                    REQUIRE(tbase.colstat[i] == LP::BStat::AtLower);
            cout << "\tChecked non-tour edges" << endl;

            AND_THEN("We can instate it in a CoreLP") {
                REQUIRE_NOTHROW(T.instate(core));
                cout << "\tInstated ActiveTour in CoreLP" << endl;

                AND_THEN("We can instate a new one from scratch") {
                    vector<int> load_tour;
                    REQUIRE_NOTHROW(util::get_tour_nodes(ncount, load_tour,
                                                         "test_data/tours/" +
                                                         prob + ".sol"));
                    double newval = 0.0;
                    vector<Graph::Edge> missing_edges;
                    for (int i = 0; i < ncount; ++i) {
                        EndPts e(load_tour[i], load_tour[(i + 1) % ncount]);
                        int len = inst.edgelen(e.end[0], e.end[1]);
                        newval += len;
                        if (graph.find_edge_ind(e.end[0], e.end[1]) == -1)
                            missing_edges.emplace_back(e.end[0], e.end[1],
                                                       len);
                    }

                    if (!missing_edges.empty())
                        REQUIRE_NOTHROW(core.add_edges(missing_edges, true));

                    cout <<"\tLoaded tour of length " << newval
                         << ", added " << missing_edges.size()
                         << " missing edges" << endl;

                    REQUIRE_NOTHROW(T = LP::ActiveTour(load_tour, core,
                                                       graph));
                    REQUIRE(T.length() == Approx(newval));
                    cout << "\tLoaded tour instated in CoreLP" << endl;
                }
            }
        }


    }
    }
}

#endif
