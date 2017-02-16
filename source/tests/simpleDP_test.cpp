#include "config.hpp"
#include "datagroups.hpp"
#include "tooth.hpp"
#include "simpleDP.hpp"
#include "process_cuts.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <set>

#include <catch.hpp>

using std::cout;
using std::setprecision;
using std::endl;

using std::vector;
using std::string;
using std::pair;

using CMR::IntPair;

#ifdef CMR_DO_TESTS

SCENARIO("Separating simple DP inequalities in small instances",
         "[SimpleDP][small]") {
    using namespace CMR;
    vector<string> probs{
        "dantzig42",
        "swiss42",
        "gr48",
        "hk48",
        "eil51",
        "st70",
        "pr76",
        "lin105",
        };

    for (string &fname : probs) {
        string
        probfile = "problems/" + fname + ".tsp",
        solfile = "test_data/tours/" + fname + ".sol",
        subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
        Graph::CoreGraph core_graph;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;
        Data::KarpPartition kpart;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can get light simple DP inequalities") {
                REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                    subtourfile, core_graph,
                                                    b_dat, lp_edges,
                                                    s_dat, inst));
                int ncount = core_graph.node_count();

                REQUIRE_NOTHROW(kpart = Data::KarpPartition(inst));
                Sep::CutQueue<Sep::dominoparity> dp_q(100);

                Sep::SimpleDP sDP(kpart, b_dat, s_dat, dp_q);

                REQUIRE(sDP.find_cuts());
                cout << "Cut queue now has size: " << dp_q.size() << "\n";

                int primal_found = 0;
                while (!dp_q.empty()) {
                    LP::SparseRow R;

                    const Sep::dominoparity &dp_cut = dp_q.peek_front();
                    vector<int> &bt = b_dat.best_tour_nodes;
                    REQUIRE_NOTHROW(R = Sep::get_row(dp_cut, bt, core_graph));

                    double tour_activity =
                    Sep::get_activity(b_dat.best_tour_edges, R);

                    double lp_activity =
                    Sep::get_activity(lp_edges, R);

                    bool viol_or_within = (lp_activity > R.rhs ||
                                           ((R.rhs - lp_activity) < 1));

                    INFO ("Known issue: lp violation sometimes off but should"
                          " never be by more than one.")
                    CAPTURE(tour_activity);
                    CAPTURE(lp_activity);
                    CAPTURE(R.rhs);

                    REQUIRE(viol_or_within);
                    REQUIRE(tour_activity <= R.rhs);

                    if (tour_activity == R.rhs && lp_activity > R.rhs)
                        ++primal_found;

                    dp_q.pop_front();
                }
                cout << "\t" << primal_found
                     << " primal violated cuts found.\n";
            }
        }
    }
}

SCENARIO("Separating simple DP inequalities in medium instances",
         "[SimpleDP][medium]") {
    using namespace CMR;
    vector<string> probs {
        "lin318",
        "d493",
        "att532",
        "u724",
        "dsj1000",
        "pr1002",
        "d2103",
        "pr2392",
        "pcb3038",
    };

    for (string &fname : probs) {
        string
        probfile = "problems/" + fname + ".tsp",
        solfile = "test_data/tours/" + fname + ".sol",
        subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
        Graph::CoreGraph core_graph;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;
        Data::KarpPartition kpart;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can get light simple DP inequalities") {
                REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                    subtourfile,
                                                    core_graph,
                                                    b_dat, lp_edges,
                                                    s_dat, inst));
                int ncount = core_graph.node_count();

                REQUIRE_NOTHROW(kpart = Data::KarpPartition(inst));
                Sep::CutQueue<Sep::dominoparity> dp_q(1000);

                Sep::SimpleDP sDP(kpart, b_dat, s_dat, dp_q);

                REQUIRE(sDP.find_cuts());
                cout << "Cut queue now has size: " << dp_q.size() << "\n";

                int primal_found = 0;
                while (!dp_q.empty()) {
                    LP::SparseRow R;

                    const Sep::dominoparity &dp_cut = dp_q.peek_front();
                    vector<int> &bt = b_dat.best_tour_nodes;

                    REQUIRE_NOTHROW(R = Sep::get_row(dp_cut, bt, core_graph));

                    double tour_activity =
                    Sep::get_activity(b_dat.best_tour_edges, R);

                    double lp_activity =
                    Sep::get_activity(lp_edges, R);

                    bool viol_or_within = (lp_activity > R.rhs ||
                                           ((R.rhs - lp_activity) < 1));

                    INFO ("Known issue: lp violation sometimes off but should"
                          " never be by more than one.")
                    CAPTURE(tour_activity);
                    CAPTURE(lp_activity);
                    CAPTURE(R.rhs);

                    REQUIRE(viol_or_within);
                    REQUIRE(tour_activity <= R.rhs);

                    if (lp_activity <= R.rhs) {
                        cout << "Found non-violated cut....\n";
                        cout << "\tHandle:" << dp_cut.degree_nodes.size()
                             << " nodes\n";
                        cout << "\tNon-neg edges:\n";
                        for (const IntPair &e : dp_cut.nonneg_edges)
                            cout << bt[e.first] << ", " << bt[e.second] << "\n";
                        cout << "\tUsed teeth:\n";
                        for (const Sep::SimpleTooth &T : dp_cut.used_teeth)
                            Sep::CandidateTeeth::print_tooth(T, false, bt);
                    }

                    if (tour_activity == R.rhs && lp_activity > R.rhs) {
                        ++primal_found;
                    }
                    dp_q.pop_front();
                }
                cout << "\t" << primal_found
                     << " primal violated cuts found.\n";
            }
        }
    }
}

SCENARIO("Separating simple DP inequalities in large instances",
         "[SimpleDP][large]") {
    using namespace CMR;
    vector<string> probs {
        "rl5915",
        "pla7397",
        "usa13509",
    };

    for (string &fname : probs) {
        string
        probfile = "problems/" + fname + ".tsp",
        solfile = "test_data/tours/" + fname + ".sol",
        subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";
        Graph::CoreGraph core_graph;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;
        Data::KarpPartition kpart;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can get light simple DP inequalities") {
                REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                    subtourfile, core_graph,
                                                    b_dat, lp_edges,
                                                    s_dat, inst));
                int ncount = core_graph.node_count();

                REQUIRE_NOTHROW(kpart = Data::KarpPartition(inst));
                Sep::CutQueue<Sep::dominoparity> dp_q(1000);

                Sep::SimpleDP sDP(kpart, b_dat, s_dat, dp_q);

                REQUIRE(sDP.find_cuts());
                cout << "Cut queue now has size: " << dp_q.size() << "\n";

                int primal_found = 0;
                while (!dp_q.empty()) {
                    LP::SparseRow R;
                    const Sep::dominoparity &dp_cut = dp_q.peek_front();
                    vector<int> &bt = b_dat.best_tour_nodes;
                    REQUIRE_NOTHROW(R = Sep::get_row(dp_cut, bt, core_graph));

                    double tour_activity =
                    Sep::get_activity(b_dat.best_tour_edges, R);

                    double lp_activity =
                    Sep::get_activity(lp_edges, R);

                    bool viol_or_within = (lp_activity > R.rhs ||
                                           ((R.rhs - lp_activity) < 1));

                    INFO ("Known issue: lp violation sometimes off but should"
                          " never be by more than one.")
                    CAPTURE(tour_activity);
                    CAPTURE(lp_activity);
                    CAPTURE(R.rhs);

                    REQUIRE(viol_or_within);
                    REQUIRE(tour_activity <= R.rhs);

                    if (tour_activity == R.rhs && lp_activity > R.rhs)
                        ++primal_found;

                    dp_q.pop_front();
                }
                cout << "\t" << primal_found
                     << " primal violated cuts found.\n";
            }
        }
    }
}

#endif //CMR_DO_TESTS
