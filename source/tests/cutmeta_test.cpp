#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "meta_sep.hpp"
#include "datagroups.hpp"
#include "solver.hpp"
#include "process_cuts.hpp"
#include "hypergraph.hpp"
#include "util.hpp"
#include "err_util.hpp"

#include <catch.hpp>

#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

SCENARIO ("Searching for cut metamorphoses",
          "[Sep][MetaCuts]") {
    using namespace CMR;
    using MetaType = Sep::MetaCuts::Type;

    vector<string> probs {
        "pr76",
        "a280",
        "lin318",
        "d493",
        "att532",
        "u724",
        "dsj1000",
        "pr1002",
        "rl1304",
        "d2103",
        "pcb3038"
        };

    for (string &prob : probs) {
    GIVEN ("A solver/cutting loop run on " + prob) {
        OutPrefs prefs;
        prefs.save_tour = false;

        Solver solver("problems/" + prob + ".tsp",
                      "test_data/tours/" + prob + ".sol",
                      prob == "a280" ? 99 : 999, prefs);

        solver.cutting_loop(false, false, true);
    THEN ("We can search for cut metamorphoses") {
        const auto &core_graph = solver.graph_info();
        const auto &active_tour = solver.active_tour();

        auto &core_lp = const_cast<LP::CoreLP &>(solver.get_core_lp());

        const auto &ext_cuts = solver.get_core_lp().external_cuts();
        const auto &lp_vec = solver.get_core_lp().lp_vec();

        vector<int> island;
        Data::SupportGroup s_dat(core_graph.get_edges(), lp_vec, island,
                                 core_graph.node_count());

        for (MetaType m : {MetaType::Decker, MetaType::Handling,
                MetaType::Teething}) {
            Sep::MetaCuts mcuts(ext_cuts, core_graph.get_edges(), active_tour,
                                s_dat);
            bool found_meta = false;
            REQUIRE_NOTHROW(found_meta = mcuts.find_cuts(m));
            if (!found_meta)
                continue;

            Sep::LPcutList &meta_q = mcuts.metacuts_q();
            const vector<int> &perm = solver.best_info().perm;
            const vector<int> &tour_edges = solver.best_info().best_tour_edges;
            for (auto it = meta_q.begin(); it; it = it->next) {
                LP::SparseRow R = Sep::get_row(*it, perm, core_graph);
                double rhs = R.rhs;
                cout << "Tour activity " << std::setprecision(2)
                     << (Sep::get_activity(tour_edges, R) - R.rhs)
                     << ", lp activity "
                     << (Sep::get_activity(lp_vec, R) - R.rhs)
                     << endl;
            }

            meta_q.clear();
            cout << "\n";
        }

        bool found_tighten = false;
        Sep::MetaCuts mcuts(ext_cuts, core_graph.get_edges(), active_tour,
                            s_dat);
        REQUIRE_NOTHROW(found_tighten = mcuts.tighten_cuts());

        if (found_tighten) {
            Sep::LPcutList &meta_q = mcuts.metacuts_q();
            const vector<int> &perm = solver.best_info().perm;
            const vector<int> &tour_edges = solver.best_info().best_tour_edges;
            for (auto it = meta_q.begin(); it; it = it->next) {
                LP::SparseRow R = Sep::get_row(*it, perm, core_graph);
                double rhs = R.rhs;
                cout << "Tour activity " << std::setprecision(2)
                     << (Sep::get_activity(tour_edges, R) - R.rhs)
                     << ", lp activity "
                     << (Sep::get_activity(lp_vec, R) - R.rhs)
                     << endl;
            }
        }


    }
    }
    }
}


#endif
