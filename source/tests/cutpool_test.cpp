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
