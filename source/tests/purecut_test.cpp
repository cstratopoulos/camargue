#include "tests.hpp"
#include "util.hpp"
#include "datagroups.hpp"
#include "core_lp.hpp"
#include "cutcontrol.hpp"
#include "cliq.hpp"
#include "timer.hpp"


#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#include <cstdlib>

#include <catch.hpp>

using std::min;
using std::max;
using std::array;
using std::vector;

using std::string;
using std::to_string;
using std::cout;

#ifdef CMR_DO_TESTS

SCENARIO ("Testing 1k random PureCut",
          "[PureCut][Solver][Random]") {
    for (int i = 0; i < 5; ++i) {
        int seed = CMR::real_zeit();

        CMR::Data::Instance inst(seed, 500, 1000000);
        CMR::Data::GraphGroup g_dat(inst);
        CMR::Data::BestGroup b_dat(inst, g_dat);
        CMR::Data::KarpPartition kpart(b_dat.perm.size(),
                                       inst.ptr(), inst.seed());

        CMR::LP::CoreLP core(g_dat, b_dat);
        CMR::TourGraph TG(b_dat.best_tour_edges,
                          g_dat.m_graph.edges,
                          b_dat.perm);

        CMR::LP::PivType piv = CMR::LP::PivType::Tour;
        bool found_cuts = true;
        int round = 0;

        CMR::Timer timer("1k" + to_string(seed));
        timer.start();

        while (true) {
            ++round;
            piv = core.primal_pivot();

            if (piv == CMR::LP::PivType::FathomedTour) {
                cout << "\t ROUND " << round << "\n";
                cout << "\t Tour fathomed optimal!\n";
                break;
            }
            if (piv == CMR::LP::PivType::Tour) {
                cout << "\t ROUND " << round << "\n";
                cout << "\t Augmented to tour of length: "
                     << core.get_objval() << "\n";
                piv = CMR::LP::PivType::Frac;
                continue;
            }
                    
            CMR::Data::SupportGroup &s_dat = core.supp_data;

            CMR::Sep::CutControl control(g_dat, b_dat, s_dat, kpart);

            if (!control.find_cuts(TG)) {
                cout << "\tROUND " << round << "\n";
                cout << "\tNo cuts found, objval: "
                     << core.get_objval() << ", stat: "
                     << CMR::LP::piv_string(piv) << "\n";
                break;
            }
                    
            core.pivot_back();

            core.add_cuts(control.seg_q);
            core.add_cuts(control.fast2m_q);
            core.add_cuts(control.blkcomb_q);
            core.add_cuts(control.connect_q);

            core.add_cuts(control.dp_q);                    
        }
        timer.stop();
        cout << "\t" << (core.num_rows() - g_dat.m_graph.node_count)
             << " cuts, " << core.num_cols() << " cols in LP.\n";
        double lp_opt = core.opt();
        double tourlen = b_dat.min_tour_value;
        double gap = 1 - (lp_opt / tourlen);
        cout << "\t Optimal obj value: " << lp_opt << ", gap: "
             << std::setprecision(10) << gap << std::setprecision(6)
             << "\n";
        timer.report(false);
        cout << "\n\n";

    }
}

SCENARIO ("Testing a prototype Purecut loop",
          "[PureCut][Solver]") {
    vector<string> probs{
//        "dantzig42", "gr48", "st70", "pr76", "rd100", "lin105"
        "a280", "lin318", "gr431", "d493", "att532", "p654", "gr666", "u724"
        //"pr1002", "d2103", "pr2392"
    };

    for (string &prob : probs) {
        GIVEN (prob) {
            THEN ("We can run a pure cut solve loop.") {
                CMR::Data::Instance inst("problems/" + prob + ".tsp",
                                         CMR::real_zeit());
                CMR::Data::GraphGroup g_dat(inst);
                CMR::Data::BestGroup b_dat(inst, g_dat);
                CMR::Data::KarpPartition kpart(b_dat.perm.size(),
                                               inst.ptr(), inst.seed());

                CMR::LP::CoreLP core(g_dat, b_dat);
                CMR::TourGraph TG(b_dat.best_tour_edges,
                                  g_dat.m_graph.edges,
                                  b_dat.perm);

                CMR::LP::PivType piv = CMR::LP::PivType::Tour;
                bool found_cuts = true;
                int round = 0;

                CMR::Timer timer(prob);
                timer.start();

                while (true) {
                    ++round;
                    piv = core.primal_pivot();

                    if (piv == CMR::LP::PivType::FathomedTour) {
                        cout << "\t ROUND " << round << "\n";
                        cout << "\t Tour fathomed optimal!\n";
                        break;
                    }
                    if (piv == CMR::LP::PivType::Tour) {
                        cout << "\t ROUND " << round << "\n";
                        cout << "\t Augmented to tour of length: "
                             << core.get_objval() << "\n";
                        piv = CMR::LP::PivType::Frac;
                        continue;
                    }
                    
                    CMR::Data::SupportGroup &s_dat = core.supp_data;

                    CMR::Sep::CutControl control(g_dat, b_dat, s_dat, kpart);

                    if (!control.find_cuts(TG)) {
                        cout << "\tROUND " << round << "\n";
                        cout << "\tNo cuts found, objval: "
                             << core.get_objval() << ", stat: "
                             << CMR::LP::piv_string(piv) << "\n";
                        break;
                    }
                    
                    core.pivot_back();

                    core.add_cuts(control.seg_q);
                    core.add_cuts(control.fast2m_q);
                    core.add_cuts(control.blkcomb_q);
                    core.add_cuts(control.connect_q);

                    core.add_cuts(control.dp_q);                    
                }
                timer.stop();
                cout << "\t" << (core.num_rows() - g_dat.m_graph.node_count)
                     << " cuts, " << core.num_cols() << " cols in LP.\n";
                cout << "\t Optimal obj value: " << core.opt() << "\n";
                timer.report(false);
                cout << "\n\n";
            }
        }
    }
}


#endif //CMR_DO_TESTS
