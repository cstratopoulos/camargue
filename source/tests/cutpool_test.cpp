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

    //trying to record seeds with lots of augmentations
    vector<ProbPair> probs{
        ProbPair("d493", 99),
        ProbPair("pr1002", 1487347814),
        ProbPair("d2103", 1487347819),
        ProbPair("pcb3038", 0),
        };

    for (ProbPair &pp : probs) {
        string &prob = pp.first;
        int seed = pp.second;
        GIVEN ("A priceless cutting loop run on " + prob) {
            THEN ("We can verify the cut prices if the cut pool is nonempty") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp", seed, prefs);
                solver.cutting_loop(false, true, true);
                Graph::CoreGraph &G =
                const_cast<Graph::CoreGraph &>(solver.graph_info());

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
                Sep::PoolCuts pool_sep(EC, edges, s_dat);

                bool found_pool = false;
                double pt = util::zeit();
                REQUIRE_NOTHROW(found_pool = pool_sep.price_cuts());
                cout << "\tPriced pool in " << (util::zeit() - pt) << "s\n";
                if (!found_pool)
                    continue;

                cout << "\tViolated cuts in pool.\n\n";
                bool found_viol = false;
                const vector<double> &slacks = pool_sep.get_slacks();
                for (int i = 0; i < slacks.size(); ++i) {
                    const Sep::HyperGraph &H = pool[i];
                    double pool_slack = slacks[i];
                    double rhs = H.get_rhs();
                    LP::SparseRow R;

                    H.get_coeffs(edges, R.rmatind, R.rmatval);

                    double activity = Sep::get_activity(s_dat.lp_vec, R);
                    double manual_slack;
                    if (H.get_sense() == 'L')
                        manual_slack = rhs - activity;
                    else if (H.get_sense() == 'G')
                        manual_slack = activity - rhs;

                    INFO ("Cut type " << H.cut_type());
                    REQUIRE(manual_slack == Approx(pool_slack));
                    if (manual_slack <= Epsilon::Cut)
                        found_viol = true;
                }
                REQUIRE(found_viol);
            }
        }
    }
}
#endif
