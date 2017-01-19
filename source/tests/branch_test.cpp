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

    using namespace CMR;

    for (string &prob : probs) {
        GIVEN ("A run of cutting loop on " + prob) {
            THEN ("If the solution is fractional we can rank edges") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   99, prefs);
                LP::PivType piv = solver.cutting_loop(false);

                if (piv == LP::PivType::Frac) {
                    LP::CoreLP &core =
                    const_cast<LP::CoreLP &>(solver.get_core_lp());
                    LP::Relaxation &rel = core;

                    vector<int> md_indices;
                    vector<double> downest;
                    vector<double> upest;

                    vector<double> x = rel.lp_vec();
                    vector<int> colstat = rel.col_stat();
                    const vector<double> &tour_edges =
                    solver.tour_basis().best_tour_edges;                    
                    double tourval = solver.best_info().min_tour_value;

                    vector<int> lw_inds;

                    for (int i = 0; i < colstat.size(); ++i)
                        if (colstat[i] == 1 && !util::var_integral(x[i]))
                            lw_inds.push_back(i);
                    
                    vector<int> sb1inds =
                        ABC::length_weighted_cands(solver.graph_info().
                                                   core_graph.get_edges(),
                                                   lw_inds, x, 5);

                    cout << "\t" << sb1inds.size()
                         << " best length weighted cands.\n";
                    for (int ind : sb1inds) {
                        cout << "Edge " << ind << ", tour "
                             << tour_edges[ind] << ", lp " << x[ind] << "\n"
                             << "\tLen "
                             << solver.graph_info().core_graph.get_edge(ind).len
                             << "\n";
                    }
                    cout << "\n";

                    const vector<int> &tourcol = solver.tour_basis().colstat;
                    const vector<int> &tourrow = solver.tour_basis().rowstat;

                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tourcol, tourrow,
                                                             sb1inds, downest,
                                                             upest, 100,
                                                             tourval));
                    vector<ABC::ScoreTuple> sb1cands;

                    REQUIRE_NOTHROW(sb1cands =
                                    ABC::ranked_cands(sb1inds, downest,
                                                      upest, 100, tourval, 2));

                    cout << "\t" << sb1cands.size() << " sb1 cands\n";
                    for (auto &st : sb1cands) {
                        int ind = st.index;
                        cout << "Edge " << ind << ", tour "
                             << tour_edges[ind] << ", lp " << x[ind] << "\n"
                             << "\tDown " << st.down_est
                             << "\tUp " << st.up_est
                             << "\tScore " << st.score << "\n\n";
                    }

                    vector<int> sb2inds{sb1cands[0].index, sb1cands[1].index};

                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tourcol, tourrow,
                                                             sb2inds, downest,
                                                             upest, 500,
                                                             tourval));

                    vector<ABC::ScoreTuple> sb2cands;
                    sb2cands = ABC::ranked_cands(sb2inds, downest, upest, 100,
                                                 tourval, 1);

                    cout << "\tWinner of sb eval\n";
                    ABC::ScoreTuple win = sb2cands[0];
                    int ind = win.index;
                    cout << "Edge " << ind << ", tour "
                         << tour_edges[ind] << ", lp " << x[ind] << "\n"
                         << "\tDown " << win.down_est
                         << "\tUp " << win.up_est 
                         << "\tScore " << win.score << "\n\n";
                        
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
