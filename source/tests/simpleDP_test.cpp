#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "datagroups.hpp"
#include "tooth.hpp"
#include "simpleDP.hpp"
#include "process_cuts.hpp"
#include "io_util.hpp"

#include <iostream>
#include <iomanip>
#include <map>
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

using TimeCuts = std::pair<double, int>;
using ProbTuple = std::tuple<string, TimeCuts, TimeCuts>;

static vector<ProbTuple> bench_probs {
    ProbTuple("pcb3038", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("fl3795", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("fnl4461", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("rl5915", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("rl5934", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("pla7397", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("rl11849", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("usa13509", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("brd14051", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("d15112", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("d18512", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    ProbTuple("pla33810", TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    // ProbTuple("pla85900",  TimeCuts(0.0, 0), TimeCuts(0.0, 0)),
    };


SCENARIO ("Benchmarking karp partitioned simple DP sep",
          "[SimpleDP][figure][table]") {
    using namespace CMR;
    for (auto &pt : bench_probs) {
        string prob = std::get<0>(pt);
        string probfile = "problems/" + prob + ".tsp";
        string solfile = "test_data/tours/" + prob + ".sol";
        string subtourfile = "test_data/subtour_lp/" + prob + ".sub.x";
    for (int i : {0, 1}) {
    GIVEN ("Data for " + prob + ", dummy part: " + std::to_string(i == 0)) {
        TimeCuts &target = (i == 0 ? std::get<1>(pt) : std::get<2>(pt));
        Graph::CoreGraph core_graph;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;
        Data::KarpPartition kpart;
        REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile, subtourfile,
                                            core_graph, b_dat, lp_edges,
                                            s_dat, inst));
        int ncount = core_graph.node_count();
        Sep::CutQueue<Sep::dominoparity> dp_q;
        LP::ActiveTour act_tour(core_graph, b_dat);
        kpart = Data::KarpPartition(inst, (i == 0), false);
        
        Sep::SimpleDP sDP(kpart, act_tour, s_dat, dp_q);
        double ft = util::zeit();
        REQUIRE(sDP.find_cuts());
        target.first = ft;

        double pfound = 0;

        while (!dp_q.empty()) {
            LP::SparseRow R;

            const Sep::dominoparity &dp_cut = dp_q.peek_front();
            vector<int> &bt = b_dat.best_tour_nodes;

            REQUIRE_NOTHROW(R = Sep::get_row(dp_cut, bt, core_graph));

            double tour_activity =
            Sep::get_activity(b_dat.best_tour_edges, R);

            REQUIRE(tour_activity <= R.rhs);

            double lp_activity =
            Sep::get_activity(lp_edges, R);

            bool viol_or_within = (lp_activity > R.rhs ||
                                   ((R.rhs - lp_activity) < 1));

            vector<int> print_elist;
            vector<double> print_ecap;
            if (tour_activity == R.rhs && lp_activity > R.rhs)
                ++pfound;

            dp_q.pop_front();
        }
        target.second = pfound;
    }
    }
    }

    THEN ("Report the data") {
        for (const auto &pt : bench_probs) {
            cout << std::get<0>(pt) << "\n";
            TimeCuts nopart = std::get<1>(pt);
            TimeCuts part = std::get<2>(pt);
            cout << "no part\t" << nopart.first << "s\t" << nopart.second
                 << " cuts\n";
            cout << "part\t" << part.first << "s\t" << part.second
                 << " cuts\n";
        }
    }
}

SCENARIO("Generating cut figures",
         "[SimpleDP][DPwitness][figure][sdp-used-edges][sdp-used-edges-part]") {
    using namespace CMR;
    vector<string> probs{
        "dsj1000",
        "pr2392",
        };

    for (string &prob : probs) {
        string probfile = "problems/" + prob + ".tsp";
        string solfile = "test_data/tours/" + prob + ".sol";
        string subtourfile = "test_data/subtour_lp/" + prob + ".sub.x";
    GIVEN ("Simple DP cuts for " + prob) {
        Graph::CoreGraph core_graph;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;
        Data::KarpPartition kpart;
        REQUIRE_NOTHROW(Data::make_cut_test(probfile, solfile,
                                            subtourfile,
                                            core_graph,
                                            b_dat, lp_edges,
                                            s_dat, inst));
        int ncount = core_graph.node_count();
        Sep::CutQueue<Sep::dominoparity> dp_q(1000);
        LP::ActiveTour act_tour(core_graph, b_dat);
        kpart = Data::KarpPartition(inst, false, true);
        Sep::SimpleDP sDP(kpart, act_tour, s_dat, dp_q);

        REQUIRE(sDP.find_cuts());

    THEN ("We can print the edge set of all cuts found") {
        int cutnum = 0;

        vector<double> used_edges(core_graph.edge_count(), 0);

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

            vector<int> print_elist;
            vector<double> print_ecap;
            if (tour_activity == R.rhs && lp_activity > R.rhs
                && R.rmatind.size() < (inst.node_count() / 2)) {
                ++cutnum;
                for (int index : R.rmatind) {
                    used_edges[index] = 1.0;
                }
            }
            dp_q.pop_front();
        }

        vector<int> print_elist;
        vector<double> print_ecap;
        for (int i = 0; i < used_edges.size(); ++i)
            if (used_edges[i] == 1.0){
                print_elist.push_back(core_graph.get_edge(i).end[0]);
                print_elist.push_back(core_graph.get_edge(i).end[1]);
                print_ecap.push_back(1.0);
            }

        cout << "\t" << cutnum << " cuts recorded" << endl;
        // util::write_xy_coords(inst.ptr()->x, inst.ptr()->y, inst.node_count(),
        //                        prob + ".xy");
        // util::write_lp_edges(print_elist, print_ecap, inst.node_count(),
        //                      prob + "-dpcuts.x");
    }
    }
    }

}

SCENARIO("Separating/illustrating tiny simple DP inequalities",
         "[SimpleDP][tiny][figure]") {
    using namespace CMR;
    string probfile = "test_data/fleisB9.tsp";
    string solfile = "test_data/tours/fleisB9.sol";
    string subtourfile = "test_data/subtour_lp/fleisB9.sub.x";
    Graph::CoreGraph core_graph;
    Data::BestGroup b_dat;
    Data::SupportGroup s_dat;
    vector<double> lp_edges;
    Data::Instance inst;
    Data::KarpPartition kpart;

    GIVEN ("The figure from Fleischer et al with a tour") {
        Data::make_cut_test(probfile, solfile, subtourfile, core_graph,
                            b_dat, lp_edges, s_dat, inst);
        kpart = Data::KarpPartition(inst);
        cout << "Best tour\n";
        for (int n : b_dat.best_tour_nodes)
            cout << n << " ";
        cout << endl;

        cout << "LP solution\n";
        for (int i = 0; i < lp_edges.size(); ++i)
            if (lp_edges[i] != 0.0)
                cout << core_graph.get_edge(i) << " " << lp_edges[i] << "\n";
        THEN ("We can get dp inequalities") {
            Sep::CutQueue<Sep::dominoparity> dp_q(100);
            LP::ActiveTour act_tour(core_graph, b_dat);

            Sep::SimpleDP sDP(kpart, act_tour, s_dat, dp_q);

            REQUIRE(sDP.find_cuts());
            cout << "\tFound " << dp_q.size() << " simple DP inequalities, "
                 << "printing..." << endl;

            while (!dp_q.empty()) {
                const auto &dp_cut = dp_q.peek_front();
                LP::SparseRow R = Sep::get_row(dp_cut,
                                               b_dat.best_tour_nodes,
                                               core_graph);
                cout << "LP activity : "
                     << Sep::get_activity(lp_edges, R) << ", tour activity: "
                     << Sep::get_activity(b_dat.best_tour_edges, R) << ", "
                     << "rhs " << R.rhs << "\n";
                cout << "Handle:\n";
                for (int n : dp_cut.degree_nodes)
                    cout << b_dat.best_tour_nodes[n] << ", ";
                cout << endl;
                cout << "Non-neg edges:\n";
                for (auto p : dp_cut.nonneg_edges)
                    cout << EndPts(b_dat.best_tour_nodes[p.first],
                                   b_dat.best_tour_nodes[p.second]) << ", ";
                cout << endl;

                for (const auto &T : dp_cut.used_teeth)
                    Sep::CandidateTeeth::print_tooth(T, true,
                                                     b_dat.best_tour_nodes);
                dp_q.pop_front();
            }

        }
    }
}


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
                LP::ActiveTour act_tour(core_graph, b_dat);
                Sep::SimpleDP sDP(kpart, act_tour, s_dat, dp_q);

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
                LP::ActiveTour act_tour(core_graph, b_dat);
                Sep::SimpleDP sDP(kpart, act_tour, s_dat, dp_q);

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
                LP::ActiveTour act_tour(core_graph, b_dat);
                Sep::SimpleDP sDP(kpart, act_tour, s_dat, dp_q);

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
