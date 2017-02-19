#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "pool_sep.hpp"
#include "datagroups.hpp"
#include "solver.hpp"
#include "process_cuts.hpp"
#include "hypergraph.hpp"
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

SCENARIO ("Experimenting with CutMonitor metrics",
          "[Sep][LP][pivot_age][tour_age][CutMonitor][experiment]") {
    using namespace CMR;
    vector<string> probs{
        "swiss42",
        "dantzig42",
        "gr48",
        "eil51",
        "pr76",
        "lin105",
        };

    for (string &prob : probs) {
    GIVEN ("Graphs/LP relaxations/Solvers for " + prob) {
        OutPrefs prefs;
        Solver solver("problems/" + prob + ".tsp", 99, prefs);


        Graph::CoreGraph &core_graph =
        const_cast<Graph::CoreGraph &>(solver.graph_info());

        Data::BestGroup &b_dat =
        const_cast<Data::BestGroup &>(solver.best_info());

        LP::CoreLP &core =
        const_cast<LP::CoreLP &>(solver.get_core_lp());

        Graph::TourGraph TG(b_dat.best_tour_edges, core_graph.get_edges(),
                            b_dat.perm);
        int ncount = core_graph.node_count();

        WHEN ("We pivot and add duplicate cuts") {
            core.primal_pivot();

            Data::SupportGroup &s_dat =
            const_cast<Data::SupportGroup &>(core.support_data());

            Sep::LPcutList bq_1;
            Sep::LPcutList bq_2;

            Sep::FastBlossoms fb_sep1(s_dat.support_elist,
                                      s_dat.support_ecap, TG, bq_1);
            Sep::FastBlossoms fb_sep2(s_dat.support_elist,
                                      s_dat.support_ecap, TG, bq_2);

            bool find1 = fb_sep1.find_cuts();
            bool find2 = fb_sep2.find_cuts();

            if (!find1 && !find2) {
                cout << "\tNo fast blossoms found" << endl;
            } else {

            REQUIRE(bq_1.size() == bq_2.size());
            int numcuts = bq_1.size();
            cout << "Found " << numcuts << " fast blossoms, adding twice"
                 << endl;

            core.pivot_back();
            core.add_cuts(bq_1);
            core.add_cuts(bq_2);
            int numrows = core.num_rows();
            vector<int> colstat;
            vector<int> rowstat;

            THEN ("We can inspect the dual variables at the next pivot") {
                core.primal_pivot();
                vector<double> piv_pi = core.pi(ncount, numrows - 1);
                cout << "Pivot pi values:\n";
                for (double pi : piv_pi)
                    cout << pi << "\n";
                cout << endl;

                core.get_base(colstat, rowstat);
                cout << "Pivot basis stats:\n";
                for (int i = ncount; i < numrows; ++i)
                    cout << rowstat[i] << "\n";
                cout << endl;
                AND_THEN ("We can inspect duals/base stats at  the tour") {
                    core.pivot_back();
                    REQUIRE(core.get_objval() == b_dat.min_tour_value);
                    vector<double> tour_pi = core.pi(ncount, numrows - 1);
                    cout << "Tour pi values:\n";
                    for (double pi : piv_pi)
                        cout << pi << "\n";
                    cout << endl;

                    core.get_base(colstat, rowstat);
                    cout << "Tour basis stats of cuts:\n";
                    for (int i = ncount; i < numrows; ++i)
                        cout << rowstat[i] << "\n";
                    cout << endl;
                }
            }
            }
        }

        THEN ("We can inspect duals after sparse cut & piv until optimal") {
            LP::PivType piv = solver.cutting_loop(false, false, true);
            if (piv != LP::PivType::FathomedTour)
                continue;
            vector<int> colstat;
            vector<int> rowstat;
            vector<int> delset(core.num_rows(), 0);

            core.get_base(colstat, rowstat);
            cout << "Optimal tour cut basic statuses:\n";
            for (int i = ncount; i < core.num_rows(); ++i) {
                cout << rowstat[i] << "\n";
                if (rowstat[i] == 1)
                    delset[i] = 1;
            }
            cout << endl;

            vector<double> tourpi = core.pi(ncount, core.num_rows() - 1);
            cout << "Optimal tour cut dual values:\n";
            for (double pi : tourpi)
                cout << pi << "\n";
            cout << endl;

            core.del_set_rows(delset);
            core.primal_opt();
            cout << "Optimal obj value after deleting basic cuts: "
                 << core.get_objval() << "\n";
            cout << endl;
        }
    }
    }
}

SCENARIO ("Pricing cuts from a cutpool",
          "[Sep][PoolCuts][price_cuts]") {
    using namespace CMR;
    using ProbPair = std::pair<string, int>;

    //trying to record seeds with lots of (or at least some) augmentations
    vector<ProbPair> probs{
        ProbPair("d493", 99),
        ProbPair("gr666", 1487433073),
        ProbPair("u724", 1487432520),
        ProbPair("pr1002", 1487347814),
        ProbPair("d2103", 1487347819),
        ProbPair("pr2392", 1487353952),
         ProbPair("pcb3038", 99),
        };

    for (ProbPair &pp : probs) {
        string &prob = pp.first;
        int seed = pp.second;
        GIVEN ("A priceless cutting loop run on " + prob) {
            THEN ("We can verify the cut prices if the cut pool is nonempty") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp", seed, prefs);

                solver.cut_sel.cutpool = false;
                solver.cutting_loop(false, true, true);

                Graph::CoreGraph &G =
                const_cast<Graph::CoreGraph &>(solver.graph_info());

                const LP::TourBasis &tbase = solver.tour_basis();

                LP::CoreLP &core =
                const_cast<LP::CoreLP &>(solver.get_core_lp());

                Sep::ExternalCuts &EC =
                const_cast<Sep::ExternalCuts &>(core.external_cuts());


                const vector<Sep::HyperGraph> &pool = EC.get_cutpool();
                if (pool.empty())
                    continue;

                cout << "\tPool has " << pool.size() << " cuts\n";

                vector<double> lp_vec = core.lp_vec();

                const vector<Graph::Edge> &edges = G.get_edges();
                vector<int> island;

                Data::SupportGroup s_dat(edges, lp_vec, island, G.node_count());
                Sep::PoolCuts pool_sep(EC, edges, tbase.best_tour_edges, s_dat);

                bool found_pool = false;
                double pt = util::zeit();
                REQUIRE_NOTHROW(found_pool = pool_sep.price_cuts());
                pt = util::zeit() - pt;
                cout << "\tPriced pool in " << pt << "s\n";

                if (found_pool)
                    cout << "\tViolated primal cuts in pool.\n\n";
                else {
                    cout << "\tNo violated primal cuts in pool.\n\n";
                }

                bool found_primal = false;
                const vector<double> &lp_slacks = pool_sep.get_lp_slacks();
                const vector<double> &tour_slacks = pool_sep.get_tour_slacks();

                for (int i = 0; i < lp_slacks.size(); ++i) {
                    const Sep::HyperGraph &H = pool[i];
                    double pool_lp_slack = lp_slacks[i];
                    double pool_tour_slack = tour_slacks[i];
                    double rhs = H.get_rhs();
                    LP::SparseRow R;

                    H.get_coeffs(edges, R.rmatind, R.rmatval);

                    double lp_activity = Sep::get_activity(s_dat.lp_vec, R);
                    double tour_activity = Sep::get_activity(tbase.
                                                             best_tour_edges,
                                                             R);
                    double manual_lp_slack;
                    double manual_tour_slack;
                    string sense_str = "\0";

                    if (H.get_sense() == 'L') {
                        manual_lp_slack = rhs - lp_activity;
                        manual_tour_slack = rhs - tour_activity;
                        sense_str = " <= ";
                    } else if (H.get_sense() == 'G') {
                        sense_str = " >= ";
                        manual_lp_slack = lp_activity - rhs;
                        manual_tour_slack = tour_activity - rhs;
                    }

                    INFO("Cut type " << H.cut_type());
                    INFO("LP: " << lp_activity << sense_str << rhs) ;
                    INFO("Tour: " << tour_activity << sense_str << rhs);
                    REQUIRE(manual_lp_slack == Approx(pool_lp_slack));
                    REQUIRE(manual_tour_slack == Approx(pool_tour_slack));
                    REQUIRE(manual_tour_slack >= 0);
                    if (manual_lp_slack <= -Epsilon::Cut &&
                        tour_activity == rhs)
                        found_primal = true;
                }
                REQUIRE(found_primal == found_pool);
            }
        }
    }
}
#endif
