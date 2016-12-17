#include "tests.hpp"
#include "core_lp.hpp"
#include "datagroups.hpp"
#include "cutcontrol.hpp"
#include "util.hpp"

#include <array>
#include <vector>
#include <string>
#include <iostream>

#include <catch.hpp>

using std::array;
using std::vector;

using std::string;
using std::cout;

#ifdef CMR_DO_TESTS

/*
SCENARIO ("Adding duplicate cuts",
          "[LP][Sep][CoreLP][CutControl][Basis]") {
    vector<string> probs{"dantzig42", "pr76"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("We add duplicate cuts and pivot back") {
                THEN ("The slack for the added cuts is in the tour basis") {
                    CMR::Data::Instance inst("problems/" + prob + ".tsp", 99);
                    CMR::Data::GraphGroup g_dat(inst);
                    CMR::Data::BestGroup b_dat(inst, g_dat);
                    CMR::LP::CoreLP core(g_dat, b_dat);

                    CMR::TourGraph TG(b_dat.best_tour_edges,
                                      g_dat.m_graph.edges, b_dat.perm);

                    REQUIRE_NOTHROW(core.primal_pivot());
                    
                    CMR::Data::SupportGroup &s_dat = core.supp_data;

                    CMR::Sep::LPcutList cutq;
                    CMR::Sep::FastBlossoms fb_sep(s_dat.support_elist,
                                                  s_dat.support_ecap, TG,
                                                  cutq);

                    REQUIRE(fb_sep.find_cuts());

                    REQUIRE_NOTHROW(core.pivot_back());

                    REQUIRE_NOTHROW(core.add_cuts(cutq));
                    REQUIRE_NOTHROW(core.add_cuts(cutq));
                    
                    int ncount = g_dat.m_graph.node_count;
                    int expect_rowcount = ncount + 4;

                    REQUIRE(core.num_rows() == expect_rowcount);

                    REQUIRE_NOTHROW(core.factor_basis());

                    vector<int> rowstat, colstat;

                    core.get_base(colstat, rowstat);
                    bool found_basis = false;
                    for (int i = ncount; i < rowstat.size(); ++i) {
                        cout << "\t Row " << i << " stat: "
                             << ((rowstat[i] == 0) ? "LOWER" : "BASIC")
                             << "\n";
                        if (rowstat[i])
                            found_basis = true;
                    }

                    REQUIRE(found_basis);
                }
            }
        }
    }
}
*/

SCENARIO ("Pivoting and adding cuts",
         "[LP][Sep][CoreLP][CutControl]") {
    vector<string> probs{"dantzig42", "st70", "pr76", "lin105",
                         "lin318", "d493", "att532", "pr1002", "rl1304",
                         "d2103", "pr2392"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("We pivot to a new solution") {
                THEN ("We can find cuts, add them, and pivot back") {
                    CMR::Data::Instance inst("problems/" + prob + ".tsp", 99);
                    CMR::Data::GraphGroup g_dat(inst);
                    CMR::Data::BestGroup b_dat(inst, g_dat);
                    CMR::LP::CoreLP core(g_dat, b_dat);
                    CMR::Data::KarpPartition kpart(b_dat.perm.size(),
                                                   inst.ptr(), 99);

                    vector<double> tour1 = core.lp_vec();

                    CMR::LP::PivType piv = core.primal_pivot();
                    cout << "Pivoted to: " << CMR::LP::piv_string(piv)
                         << ", piv val: " << core.get_objval() << "\n";

                    vector<double> piv1 = core.lp_vec();

                    CMR::Data::SupportGroup &s_dat = core.supp_data;

                    CMR::TourGraph TG(b_dat.best_tour_edges,
                                      g_dat.m_graph.edges,
                                      b_dat.perm);

                    CMR::Sep::CutControl control(g_dat, b_dat, s_dat, kpart);

                    REQUIRE(control.find_cuts(TG));

                    REQUIRE_NOTHROW(core.pivot_back());

                    vector<double> tour2 = core.lp_vec();

                    REQUIRE(tour1 == tour2);


                    std::array<CMR::Sep::LPcutList*,
                               4> qlist{&control.seg_q, &control.fast2m_q,
                                        &control.blkcomb_q,
                                        &control.connect_q};

                    for (auto qptr : qlist)
                        REQUIRE_NOTHROW(core.add_cuts(*qptr));

                    REQUIRE_NOTHROW(core.add_cuts(control.dp_q));

                    piv = core.primal_pivot();
                    cout << "Pivoted to: " << CMR::LP::piv_string(piv)
                         << ", piv val: " << core.get_objval() << "\n";

                    vector<double> piv2 = core.lp_vec();

                    REQUIRE(piv1 != piv2);

                    
                    cout << "\n\n";
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
