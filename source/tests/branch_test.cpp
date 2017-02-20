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

extern "C" {
#include <concorde/INCLUDE/linkern.h>
}

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
        //"pr76",
        "a280",
        "lin318",
        "p654",
        };

    for (string &prob : probs) {
        GIVEN ("A Solver for " + prob) {
            THEN ("We can run a contra fixing ABC") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp", 99, prefs);
                LP::PivType piv = LP::PivType::Frac;

                REQUIRE_NOTHROW(piv = solver.abc(true));
                cout << "\n\tTerminated abc search: " << piv << "\n";
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
                int ncount = solver.inst_info().node_count();
                LP::PivType piv = solver.cutting_loop(ncount < 100, true, true);

                if (piv == LP::PivType::Frac) {
                    LP::CoreLP &core =
                    const_cast<LP::CoreLP &>(solver.get_core_lp());
                    LP::Relaxation &rel = core;

                    std::unique_ptr<ABC::Brancher> branch;

                    const vector<Graph::Edge> &edges = solver.graph_info().
                    get_edges();

                    const LP::TourBasis &tbase = solver.tour_basis();

                    const double tourlen = solver.best_info().min_tour_value;
                    vector<double> frac_vec = rel.lp_vec();

                    REQUIRE_NOTHROW(branch =
                                    util::make_unique<ABC::Brancher>(rel,
                                                                     edges,
                                                                     tbase,
                                                                     tourlen,
                                                                     ABC::ContraStrat::Fix));

                    auto obj = branch->next_branch_obj();
                    int ind = obj.index;
                    cout << "\tIndex " << ind << " should be next branch.\n";

                    double tour_entry = tbase.best_tour_edges[ind];
                    double lp_entry = frac_vec[ind];

                    cout << "\tLP " << lp_entry << "\tTour " << tour_entry
                         << "\n";


                    core.pivot_back(false);
                }
            }
        }
    }
}

SCENARIO ("Computing branching edges",
          "[ABC][primal_strong_branch][ranked_cands][best_estimate]") {
    vector<string> probs{
        "dantzig42",
        "pr76",
        "pr152",
        "a280",
        "lin318",
        "d493",
        "att532",
        "pr1002",
        "rl1304",
        "d2103",
        "pr2392",
        "pcb3038",
    };

    using namespace CMR;

    for (string &prob : probs) {
        GIVEN ("A run of cutting loop on " + prob) {
            THEN ("If the solution is fractional we can rank edges") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   99, prefs);
                int ncount = solver.inst_info().node_count();
                LP::PivType piv = solver.cutting_loop(ncount < 100, false,
                                                      true);

                if (piv == LP::PivType::Frac) {
                    LP::CoreLP &core =
                    const_cast<LP::CoreLP &>(solver.get_core_lp());
                    LP::Relaxation &rel = core;

                    int avgit = core.avg_itcount();

                    cout << "\nAverage number of pivots per nd pivot: "
                         << avgit << "\n\n";

                    vector<int> md_indices;
                    vector<ABC::ScorePair> downest;
                    vector<ABC::ScorePair> upest;

                    vector<double> x = rel.lp_vec();
                    vector<int> colstat = rel.col_stat();

                    vector<double> tour_edges = solver.tour_basis()
                    .best_tour_edges;
                    double tourval = solver.best_info().min_tour_value;

                    vector<int> lw_inds;

                    for (int i = 0; i < colstat.size(); ++i)
                        if (colstat[i] == 1 && !util::var_integral(x[i]))
                            lw_inds.push_back(i);

                    vector<int> sb1inds =
                        ABC::length_weighted_cands(solver.graph_info().
                                                   get_edges(),
                                                   lw_inds, x, 5);

                    // cout << "\t" << sb1inds.size()
                    //      << " best length weighted cands.\n";
                    // for (int ind : sb1inds) {
                    //     cout << "Edge " << ind << ", tour "
                    //          << tour_edges[ind] << ", lp " << x[ind] << "\n"
                    //          << "\tLen "
                    //          << solver.graph_info().core_graph.get_edge(ind).len
                    //          << "\n";
                    // }
                    // cout << "\n";

                    const vector<int> &tourcol = solver.tour_basis().colstat;
                    const vector<int> &tourrow = solver.tour_basis().rowstat;
                    vector<LP::Basis> cbases;

                    int sb1ic = std::max(100, 2 * avgit);
                    cout << "Calling primal SB with ic " << sb1ic << "\n";
                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tourcol, tourrow,
                                                             sb1inds, downest,
                                                             upest,
                                                             cbases,
                                                             sb1ic,
                                                             tourval));
                    vector<ABC::ScoreTuple> sb1cands;

                    REQUIRE_NOTHROW(sb1cands =
                                    ABC::ranked_cands(sb1inds, downest,
                                                      upest, cbases, 100,
                                                      tourval, 2));

                    // cout << "\t" << sb1cands.size() << " sb1 cands\n";
                    // for (auto &st : sb1cands) {
                    //     int ind = st.index;
                    //     cout << "Edge " << ind << ", tour "
                    //          << tour_edges[ind] << ", lp " << x[ind] << ", "
                    //          << "Priority " << st.score_priority << "\n"
                    //          << "\tDown " << st.down_est.first << " -- "
                    //          << st.down_est.second << "\n"
                    //          << "\tUp " << st.up_est.first << " -- "
                    //          << st.up_est.second << "\n"
                    //          << "\tScore " << st.score << "\n\n";
                    // }

                    vector<LP::Basis> sb2bases;
                    vector<int> sb2inds;

                    for (auto &t : sb1cands) {
                        sb2inds.push_back(t.index);
                        sb2bases.emplace_back(std::move(t.contra_base));
                    }

                    int sb2ic = std::min(5 * avgit, 500);
                    cout << "Calling primal SB2 with ic " << sb2ic << "\n";
                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tourcol, tourrow,
                                                             sb2inds, downest,
                                                             upest, sb2bases,
                                                             sb2ic,
                                                             tourval));

                    vector<ABC::ScoreTuple> sb2cands;
                    sb2cands = ABC::ranked_cands(sb2inds, downest, upest,
                                                 sb2bases, 100, tourval, 1);

                    cout << "\tWinner of sb eval\n";
                    ABC::ScoreTuple &win = sb2cands[0];
                    int ind = win.index;
                    cout << "Edge " << ind << ", tour "
                         << tour_edges[ind] << ", lp " << x[ind] << ", "
                         << "Priority " << win.score_priority << "\n"
                         << "\tDown " << win.down_est.first << " -- "
                         << win.down_est.second << "\n"
                         << "\tUp " << win.up_est.first << " -- "
                         << win.up_est.second << "\n"
                         << "\tScore " << win.score << "\n\n";

                    vector<int> elist;
                    vector<int> elen;

                    solver.graph_info().get_elist(elist, elen);

                    Data::Instance sp_inst(prob, 99, ncount, elist, elen);

                    int default_len = sp_inst.ptr()->default_len;

                    if (tour_edges[ind] == 1)
                        elen[ind] = default_len;
                    else
                        elen[ind] = - default_len;

                    sp_inst = Data::Instance(prob, 99, ncount, elist, elen);

                    vector<int> tour_nodes = solver.best_info()
                    .best_tour_nodes;
                    vector<int> out_cyc(ncount);

                    CCrandstate rstate;
                    CCutil_sprand(99, &rstate);
                    double val;

                    cout << "\n\tRunning branch linkern....\n";

                    REQUIRE_FALSE(CClinkern_tour(ncount, sp_inst.ptr(),
                                                 elen.size(), &elist[0],
                                                 ncount, 500,
                                                 // &tour_nodes[0],
                                                 NULL,
                                                 &out_cyc[0],
                                                 &val,
                                                 1, 0, 0,
                                                 (char *) NULL,
                                                 CC_LK_GEOMETRIC_KICK,
                                                 &rstate));

                    const auto &graph = solver.graph_info();
                    const auto &edges = graph.get_edges();


                    bool found_branch_edge = false;
                    double newval = 0.0;

                    vector<int> new_tour_edges(tour_edges.size(), 0);

                    for (int i = 0; i < ncount; ++i) {
                        EndPts e(out_cyc[i], out_cyc[(i + 1) % ncount]);
                        newval += solver.inst_info().edgelen(e.end[0],
                                                             e.end[1]);
                        int ii = graph.find_edge_ind(e.end[0], e.end[1]);
                        if (ii == -1) {
                            tour_edges.push_back(0);
                            new_tour_edges.push_back(1);
                            cout << "Sparse lk added an edge??\n";
                        } else
                            new_tour_edges[ii] = 1;
                    }

                    REQUIRE(tour_edges[ind] != new_tour_edges[ind]);

                    cout << "\tRaw lk val: " << val << "\n";

                    double orig_val = solver.best_info().min_tour_value;
                    cout << "\tOrig val: " << orig_val
                         << "\n\tbranch tour: " << newval << "\n\n";

                    if (newval < orig_val)
                        cout << "\t||||| BETTER TOUR!!! |||||\n\n";

                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
