#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "lp_interface.hpp"
#include "branch_util.hpp"
#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

#include <cstdlib>

#include <catch.hpp>

using std::min;
using std::max;
using std::abs;
using std::array;
using std::vector;
using std::pair;

using std::string;
using std::to_string;
using std::cout;

SCENARIO ("Computing branching edges",
          "[ABC]") {
    vector<string> probs{
        "pr152",
        "a280",
        "lin318",
    };

    for (string &prob : probs) {
        GIVEN ("A run of cutting loop on " + prob) {
            THEN ("If the solution is fractional we can rank edges") {
                CMR::OutPrefs prefs;
                CMR::Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   99, prefs);
                CMR::LP::PivType piv = solver.cutting_loop(false);

                if (piv == CMR::LP::PivType::Frac) {
                    CMR::LP::CoreLP &core =
                    const_cast<CMR::LP::CoreLP &>(solver.get_core_lp());
                    CMR::LP::Relaxation &rel = core;

                    vector<int> colstat;
                    vector<int> rowstat;
                    vector<double> x = rel.lp_vec();
                    vector<int> md_indices;
                    
                    rel.get_base(colstat, rowstat);

                    for (int i = 0; i < colstat.size(); ++i)
                        if (colstat[i] == CMR::LP::BStat::Basic &&
                            !CMR::util::var_integral(x[i]))
                            md_indices.push_back(i);

                    REQUIRE_FALSE(md_indices.empty());

                    vector<double> downest;
                    vector<double> upest;

                    REQUIRE_NOTHROW(rel.get_penalties(md_indices, downest,
                                                      upest));
                    
                    vector<CMR::ABC::ScorePair> md_cands;

                    REQUIRE_NOTHROW(md_cands =
                                    CMR::ABC::ranked_cands(md_indices, downest,
                                                           upest, 10, 5));
                    const vector<double> &tour_edges =
                    solver.tour_basis().best_tour_edges;
                    
                    double tourval = solver.best_info().min_tour_value;

                    vector<int> sb1inds;

                    cout << "\t" << md_cands.size() << " driebek cands\n";
                    
                    for (auto &sp : md_cands) {
                        int ind = sp.first;
                        sb1inds.push_back(ind);
                        cout << "\tEdge " << ind << "\n\t\ttour: "
                             << tour_edges[ind] << "\t"
                             << "lp: " << x[ind] << "\t"
                             << "Est: " << sp.second << "\n";
                    }
                    cout << "\n";

                    const vector<int> &tourcol = solver.tour_basis().colstat;
                    const vector<int> &tourrow = solver.tour_basis().rowstat;

                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tourcol, tourrow,
                                                             sb1inds, downest,
                                                             upest, 100,
                                                             tourval));
                    vector<CMR::ABC::ScorePair> sb1cands;

                    REQUIRE_NOTHROW(sb1cands =
                                    CMR::ABC::ranked_cands(sb1inds, downest,
                                                           upest, 100, 2));

                    cout << "\t" << sb1cands.size() << " sb1 cands\n";
                    for (auto &sp : sb1cands) {
                        int ind = sp.first;
                        sb1inds.push_back(ind);
                        cout << "\tEdge " << ind << "\n\t\ttour: "
                             << tour_edges[ind] << "\t"
                             << "lp: " << x[ind] << "\t"
                             << "Est: " << sp.second << "\n";
                    }
                    cout << "\n";                    
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
