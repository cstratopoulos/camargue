#include "tests.hpp"
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

constexpr double EpsZero = CMR::Epsilon::Zero;
using ScorePair = std::pair<int, double>;

static bool sp_sort_greater(const ScorePair &p, const ScorePair &r)
{
    return p.second > r .second;
}

static bool is_integral(double d) { return d >= EpsZero && d <= 1 - EpsZero; }

#ifdef CMR_DO_TESTS

SCENARIO ("Computing branching edges",
          "[Solver][Branch][ABC]") {
    vector<string> probs{
        "d493", "pr1002", "d2103", "pr2392"
    };

    for (string &prob : probs) {
        GIVEN ("A run of cutting loop on " + prob) {
            THEN ("If the solution is fractional we can rank edges") {
                CMR::OutPrefs prefs;
                CMR::Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   99, prefs);

                CMR::LP::PivType piv;

                REQUIRE_NOTHROW(piv = solver.cutting_loop());

                if (piv == CMR::LP::PivType::Frac) {
                    vector<int> frac_colstat;
                    vector<int> frac_rowstat;
                    vector<double> frac_x;

                    CMR::LP::Relaxation &rel = solver.relax();
                    
                    const CMR::LP::TourBasis &tour_base = solver.tour_basis();

                    vector<double> tour_edges = tour_base.best_tour_edges;
                    vector<int> tour_colstat = tour_base.colstat;
                    vector<int> tour_rowstat = tour_base.rowstat;
                    
                    REQUIRE_NOTHROW(rel.get_base(frac_colstat, frac_rowstat));
                    REQUIRE_NOTHROW(rel.get_x(frac_x));

                    REQUIRE(frac_colstat.size() > 0);

                    vector<int> cand_inds;
                    vector<double> downratio;
                    vector<double> upratio;

                    for (int i = 0; i < frac_x.size(); ++i)
                        if(is_integral(frac_x[i]) &&
                            frac_colstat[i] == CPX_BASIC)
                            cand_inds.push_back(i);

                    REQUIRE(cand_inds.size() > 0);
                    REQUIRE_NOTHROW(rel.get_penalties(cand_inds, downratio,
                                                      upratio));
                    
                    vector<ScorePair> db_scores;

                    for (int i = 0; i < cand_inds.size(); ++i) {
                        ScorePair p(cand_inds[i],
                                    CMR::LP::branch_score(10, downratio[i],
                                                          upratio[i]));
                        db_scores.push_back(p);
                    }

                    std::sort(db_scores.begin(), db_scores.end(),
                              sp_sort_greater);

                    db_scores.resize(5);
                    for (ScorePair &p : db_scores)
                        cout << "Var " << p.first << " has score "
                             << p.second << "\n";

                    cand_inds.clear();

                    for (ScorePair &p : db_scores)
                        cand_inds.push_back(p.first);

                    vector<double> downobj;
                    vector<double> upobj;

                    const CMR::Data::BestGroup &b_dat = solver.best_info();

                    double tourval = b_dat.min_tour_value;

                    cout << "Doing 100 itlim primal sb with 5 cands...";
                    double sb1 = CMR::zeit();
                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tour_colstat,
                                                             tour_rowstat,
                                                             cand_inds,
                                                             downobj,
                                                             upobj, 100,
                                                             tourval));
                    cout << "Done in " << (CMR::zeit() - sb1) << "s\n";

                    vector<double> after_x;
                    
                    REQUIRE_NOTHROW(after_x = rel.lp_vec());
                    REQUIRE(after_x == tour_edges);
                    

                    vector<ScorePair> sb1_scores;

                    for (int i = 0; i < cand_inds.size(); ++i) {
                        ScorePair p(cand_inds[i],
                                    CMR::LP::branch_score(100, downobj[i],
                                                          upobj[i]));
                        sb1_scores.push_back(p);
                    }

                    std::sort(sb1_scores.begin(), sb1_scores.end(),
                              sp_sort_greater);

                    for (ScorePair &p : sb1_scores)
                        cout << "Var " << p.first << " has sb1 score "
                             << p.second << "\n";

                    vector<int> sb1_5_cands{sb1_scores[0].first,
                                            sb1_scores[1].first};

                    double sb2 = CMR::zeit();

                    cout << "500 itlim primal sb w 2 cands...";
                    REQUIRE_NOTHROW(rel.primal_strong_branch(tour_edges,
                                                             tour_colstat,
                                                             tour_rowstat,
                                                             sb1_5_cands,
                                                             downobj,
                                                             upobj, 500,
                                                             tourval));
                    cout << "Done in " << (CMR::zeit() - sb2) << "s\n";

                    vector<ScorePair> sb2_scores;
                    
                    for (int i = 0; i < 2; ++i) {
                        ScorePair p(sb1_5_cands[i],
                                  CMR::LP::branch_score(100, downobj[i],
                                                        upobj[i]));
                        sb2_scores.push_back(p);
                    }

                    for (ScorePair &p : sb2_scores)
                        cout << "Var " << p.first << " has sb2 score "
                             << p.second << "\n";

                    rel.copy_start(frac_x, frac_colstat, frac_rowstat);
                    rel.factor_basis();

                    cout << "Doing CPLEX dual sb 100 itlim 5 cands NO opt....";
                    double dsb = CMR::zeit();

                    REQUIRE_NOTHROW(rel.dual_strong_branch(cand_inds,
                                                           downobj, upobj,
                                                           100, tourval));
                    cout << "Done in " << (CMR::zeit() - dsb) << "s\n";

                    vector<ScorePair> dsb_scores;

                    for (int i = 0; i < cand_inds.size(); ++i) {
                        ScorePair p(cand_inds[i],
                                    CMR::LP::branch_score(100, downobj[i],
                                                          upobj[i]));
                        dsb_scores.push_back(p);
                    }

                    cout << "\tRaw dual sb data\n";
                    for (int i = 0; i < cand_inds.size(); ++i) {
                        cout << "var " << cand_inds[i] << " downobj: "
                             << downobj[i] << ", upobj: "
                             << upobj[i] << "tour vec: "
                             << tour_edges[i] << "\n";
                    }

                    std::sort(dsb_scores.begin(), dsb_scores.end(),
                              sp_sort_greater);

                    for (ScorePair &p : dsb_scores)
                        cout << "Var " << p.first << " has sb2 score "
                             << p.second << "\n";
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
