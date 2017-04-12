#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "separator.hpp"
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

SCENARIO ("Adding duplicate cuts to the pool",
          "[Sep][ExternalCuts][pool_add]") {
    using namespace CMR;
    vector<string> probs{
        "pr76",
        "d493",
        "rl1304",
        };

    for (string &prob : probs) {
    GIVEN ("Subtour LP data for " + prob) {
        string
        probfile = "problems/" + prob + ".tsp",
        solfile = "test_data/tours/" + prob + ".sol",
        subtourfile = "test_data/subtour_lp/" + prob + ".sub.x";

        Graph::CoreGraph core_graph;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        std::vector<double> lp_edges;
        Sep::LPcutList cutq;
        Data::make_cut_test(probfile, solfile, subtourfile, core_graph, b_dat,
                            lp_edges, s_dat);

    WHEN ("Cuts are found") {
        vector<double> d_tour_edges(b_dat.best_tour_edges.begin(),
                                    b_dat.best_tour_edges.end());
        Sep::TourGraph TG(d_tour_edges,
                          core_graph.get_edges(),
                          b_dat.perm);
        for (int &i : s_dat.support_elist) i = b_dat.perm[i];

        Sep::FastBlossoms fb_sep(s_dat.support_elist,
                                 s_dat.support_ecap, TG, cutq);
        if (!fb_sep.find_cuts())
            continue;

        THEN ("Adding the same cut twice doesn't increase pool size") {
            CCtsp_lpcuts *cc_pool;
            auto pguard = util::make_guard([&cc_pool]()
                                           { CCtsp_free_cutpool(&cc_pool); });

            int ncount = core_graph.node_count();
            REQUIRE_FALSE(CCtsp_init_cutpool(&ncount, NULL, &cc_pool));
            REQUIRE(cc_pool->cutcount == 0);

            for (CCtsp_lpcut_in *c = cutq.begin(); c; c = c->next)
                REQUIRE_FALSE(CCtsp_add_to_cutpool_lpcut_in(cc_pool, c));


            REQUIRE(cc_pool->cutcount == cutq.size());

            for (CCtsp_lpcut_in *c = cutq.begin(); c; c = c->next)
                REQUIRE_FALSE(CCtsp_add_to_cutpool_lpcut_in(cc_pool, c));


            REQUIRE(cc_pool->cutcount == cutq.size());

            AND_THEN ("Re-separating from pool also doesn't increase size") {
                Sep::LPcutList poolq;
                Sep::PoolCuts pool_sep(s_dat.support_elist, s_dat.support_ecap,
                                       TG, poolq, cc_pool, 99);
                if (!pool_sep.find_cuts())
                    continue;

                for (CCtsp_lpcut_in *c = poolq.begin(); c; c = c->next)
                    REQUIRE_FALSE(CCtsp_add_to_cutpool_lpcut_in(cc_pool, c));

                REQUIRE(cc_pool->cutcount == cutq.size());
            }

        }
    }
    }
    }
}

SCENARIO ("Experimenting with CutMonitor metrics",
          "[Sep][LP][pivot_age][CutMonitor][experiment]") {
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

        Sep::TourGraph TG(solver.active_tour().edges(), core_graph.get_edges(),
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

            core.pivot_back(false);
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
                    core.pivot_back(false);
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
          "[Sep][PoolCuts][pool_sep][tighten_pool]") {
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

                LP::CoreLP &core =
                const_cast<LP::CoreLP &>(solver.get_core_lp());

                Sep::ExternalCuts &EC =
                const_cast<Sep::ExternalCuts &>(core.external_cuts());

                core.primal_opt();

                vector<double> lp_vec = core.lp_vec();

                const vector<Graph::Edge> &edges = G.get_edges();
                vector<int> island;

                Data::SupportGroup s_dat(edges, lp_vec, island, G.node_count());
                Data::KarpPartition kpart;
                Sep::Separator sep(G.get_edges(), solver.active_tour(),
                                   s_dat, kpart, 99);
                sep.verbose = true;

                bool found_pool = false;
                REQUIRE_NOTHROW(found_pool = sep.pool_sep(EC));

                if (found_pool)
                    cout << "\tViolated primal cuts in pool.\n\n";
                else {
                    cout << "\tNo violated primal cuts in pool.\n\n";
                }

            AND_THEN("We can search for tighten pool cuts") {
                bool found_tight = false;
                REQUIRE_NOTHROW(found_tight = sep.tighten_pool(EC));
                AND_THEN("We can search for consec1s combs") {
                    bool found_con1 = false;
                    REQUIRE_NOTHROW(found_con1 = sep.consec1_sep(EC));
                }
            }
            }
        }
    }
}
#endif
