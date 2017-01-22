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
    vector<string> probs {"dantzig42", "swiss42", "gr48", "hk48", "eil51",
                          "st70", "pr76", "lin105"};

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
        Data::KarpPartition kpart;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can get light simple DP inequalities") {
                REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                         subtourfile, g_dat,
                                                         b_dat, lp_edges,
                                                         s_dat, inst));
                int ncount = g_dat.core_graph.node_count();

                REQUIRE_NOTHROW(kpart = Data::KarpPartition(ncount,
                                                                 inst.ptr(),
                                                                 99));
                Sep::CutTranslate translator(g_dat);
                Sep::CutQueue<Sep::dominoparity> dp_q(100);

                Sep::SimpleDP sDP(g_dat, kpart, b_dat, s_dat, dp_q);

                REQUIRE(sDP.find_cuts());
                cout << "Cut queue now has size: " << dp_q.size() << "\n";

                int primal_found = 0;
                while (!dp_q.empty()) {
                    vector<int> rmatind;
                    vector<double> rmatval;
                    char sense;
                    double rhs;
	  
                    const Sep::dominoparity &dp_cut = dp_q.peek_front();
                    vector<int> &bt = b_dat.best_tour_nodes;
                    double tour_activity, lp_activity;
                    REQUIRE_NOTHROW(translator.get_sparse_row(dp_cut, bt,
                                                              rmatind,
                                                              rmatval, sense,
                                                              rhs));

                    translator.get_activity(tour_activity,
                                            b_dat.best_tour_edges,
                                            rmatind, rmatval);
                    translator.get_activity(lp_activity, lp_edges, rmatind,
                                            rmatval);

                    CHECK(lp_activity > rhs);
                    REQUIRE(tour_activity <= rhs);
	  
                    if (tour_activity == rhs && lp_activity > rhs)
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
        "d493", "att532", "u724",
        "dsj1000", "pr1002",
        "d2103", "pr2392",
        "pcb3038",
    };

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
        Data::KarpPartition kpart;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can get light simple DP inequalities") {
                REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                         subtourfile, g_dat,
                                                         b_dat, lp_edges,
                                                         s_dat, inst));
                int ncount = g_dat.core_graph.node_count();

                REQUIRE_NOTHROW(kpart = Data::KarpPartition(ncount,
                                                                 inst.ptr(),
                                                                 99));
                Sep::CutTranslate translator(g_dat);
                Sep::CutQueue<Sep::dominoparity> dp_q(1000);

                Sep::SimpleDP sDP(g_dat, kpart, b_dat, s_dat, dp_q);

                REQUIRE(sDP.find_cuts());
                cout << "Cut queue now has size: " << dp_q.size() << "\n";

                int primal_found = 0;
                while (!dp_q.empty()) {
                    vector<int> rmatind;
                    vector<double> rmatval;
                    char sense;
                    double rhs;
	  
                    const Sep::dominoparity &dp_cut = dp_q.peek_front();
                    vector<int> &bt = b_dat.best_tour_nodes;
                    double tour_activity, lp_activity;
          
                    REQUIRE_NOTHROW(translator.get_sparse_row(dp_cut, bt,
                                                              rmatind,
                                                              rmatval, sense,
                                                              rhs));

                    translator.get_activity(tour_activity,
                                            b_dat.best_tour_edges,
                                            rmatind, rmatval);
                    translator.get_activity(lp_activity, lp_edges, rmatind,
                                            rmatval);

                    CHECK(lp_activity > rhs);
                    REQUIRE(tour_activity <= rhs);

                    if (lp_activity <= rhs) {
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
	  
                    if (tour_activity == rhs && lp_activity > rhs) {
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
        "rl5915", "pla7397",
        "usa13509"
    };

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
        Data::KarpPartition kpart;

        GIVEN("A subtour polytope LP solution for " + fname) {
            THEN("We can get light simple DP inequalities") {
                REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                                         subtourfile, g_dat,
                                                         b_dat, lp_edges,
                                                         s_dat, inst));
                int ncount = g_dat.core_graph.node_count();

                REQUIRE_NOTHROW(kpart = Data::KarpPartition(ncount,
                                                                 inst.ptr(),
                                                                 99));
                Sep::CutTranslate translator(g_dat);
                Sep::CutQueue<Sep::dominoparity> dp_q(1000);

                Sep::SimpleDP sDP(g_dat, kpart, b_dat, s_dat, dp_q);

                REQUIRE(sDP.find_cuts());
                cout << "Cut queue now has size: " << dp_q.size() << "\n";

                int primal_found = 0;
                while (!dp_q.empty()) {
                    vector<int> rmatind;
                    vector<double> rmatval;
                    char sense;
                    double rhs;
	  
                    const Sep::dominoparity &dp_cut = dp_q.peek_front();
                    vector<int> &bt = b_dat.best_tour_nodes;
                    double tour_activity, lp_activity;
                    REQUIRE_NOTHROW(translator.get_sparse_row(dp_cut, bt,
                                                              rmatind,
                                                              rmatval, sense,
                                                              rhs));

                    translator.get_activity(tour_activity,
                                            b_dat.best_tour_edges,
                                            rmatind, rmatval);
                    translator.get_activity(lp_activity, lp_edges, rmatind,
                                            rmatval);

                    CHECK(lp_activity > rhs);
                    REQUIRE(tour_activity <= rhs);
	  
                    if (tour_activity == rhs && lp_activity > rhs)
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
