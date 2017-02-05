#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "lp_interface.hpp"
#include "brancher.hpp"
#include "branch_util.hpp"
#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <memory>
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

SCENARIO ("Running a Solver with contra Fix Brancher",
          "[ABC][Brancher][ContraStrat][Fix]") {
    using namespace CMR;
    vector<string> probs{
        "dantzig42",
        "pr76",
        "a280",
        "lin318"
        };

    for (string &prob : probs) {
        GIVEN ("A Solver for " + prob) {
            THEN ("We can run a contra fixing ABC") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp", 99, prefs);
                LP::PivType piv = LP::PivType::Frac;

                REQUIRE_NOTHROW(piv = solver.abc(true));
                cout << piv << "\n";
            }
        }
    }
    
}

SCENARIO ("Instating a Brancher and getting problems",
          "[ABC][Brancher]") {
    vector<string> probs{
        "dantzig42",
        "pr152",
        "a280",
        "lin318",
        "d493",
        };

    using namespace CMR;

    for (string &prob : probs) {
        GIVEN ("A cutting loop run on " + prob) {
            THEN ("We can get a first prob from Brancher and try to enforce it")
            {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp",
                              //prob + ".sol",
                              99, prefs);
                LP::PivType piv = solver.cutting_loop(solver.inst_info().node_count() < 100, true);

                if (piv == LP::PivType::Frac) {
                    LP::CoreLP &core =
                    const_cast<LP::CoreLP &>(solver.get_core_lp());
                    LP::Relaxation &rel = core;

                    std::unique_ptr<ABC::Brancher> branch;
                    
                    const vector<Graph::Edge> &edges = solver.graph_info().
                    core_graph.get_edges();

                    const LP::TourBasis &tbase = solver.tour_basis();

                    const double tourlen = solver.best_info().min_tour_value;
                    vector<double> frac_vec = rel.lp_vec();

                    REQUIRE_NOTHROW(branch =
                                    util::make_unique<ABC::Brancher>(rel,
                                                                     edges,
                                                                     tbase,
                                                                     tourlen,
                                                                     ABC::ContraStrat::Fix));

                    int ind = branch->branch_edge_index();
                    cout << "\tIndex " << ind << " should be next branch.\n";
                    
                    double tour_entry = tbase.best_tour_edges[ind];
                    double lp_entry = frac_vec[ind];

                    cout << "\tLP " << lp_entry << "\tTour " << tour_entry
                         << "\n";

                    
                    double zero_coeff;
                    double one_coeff;

                    ABC::dive_coeffs(tourlen, zero_coeff, one_coeff);
                    cout << "\tSuggested dive coeffs\t" << zero_coeff << "\t"
                         << one_coeff << "\n";

                    core.pivot_back();
                    double changeto = tour_entry == 0 ? one_coeff : zero_coeff;
                    
                    REQUIRE_NOTHROW(rel.change_obj(ind, changeto));
                    core.primal_pivot();
                    double new_lp_entry = rel.lp_vec()[ind];
                    REQUIRE(new_lp_entry != tour_entry);
                    cout << "\tDived objval: " << rel.get_objval() << "\n\n";
                }
            }
        }
    }
}

SCENARIO ("Computing branching edges",
          "[.ABC]") {
    vector<string> probs{
        "pr152",
        "a280",
        "lin318",
        "d493",
        "pr1002",
    };

    using namespace CMR;

    for (string &prob : probs) {
        GIVEN ("A run of cutting loop on " + prob) {
            THEN ("If the solution is fractional we can rank edges") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   99, prefs);
                LP::PivType piv = solver.cutting_loop(false, true);

                if (piv == LP::PivType::Frac) {
                    LP::CoreLP &core =
                    const_cast<LP::CoreLP &>(solver.get_core_lp());
                    LP::Relaxation &rel = core;

                    vector<int> md_indices;
                    vector<ABC::ScorePair> downest;
                    vector<ABC::ScorePair> upest;

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
                             << tour_edges[ind] << ", lp " << x[ind] << ", "
                             << "Priority " << st.score_priority << "\n"
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
                         << tour_edges[ind] << ", lp " << x[ind] << ", "
                         << "Priority " << win.score_priority << "\n"
                         << "\tDown " << win.down_est
                         << "\tUp " << win.up_est 
                         << "\tScore " << win.score << "\n\n";
                        
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
