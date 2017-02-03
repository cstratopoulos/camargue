#include "config.hpp"

#ifdef CMR_DO_TESTS
#if CMR_HAVE_SAFEGMI

#include "solver.hpp"
#include "safeGMI.hpp"
#include "util.hpp"
#include "timer.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>




#include <cstdlib>

#include <catch.hpp>

using std::min;
using std::max;

using std::array;
using std::vector;
using std::unique_ptr;

using std::string;
using std::to_string;
using std::cout;



SCENARIO ("Generating safe Gomory cuts",
          "[Sep][SafeGomory]") {
    using namespace CMR;
    vector<string> probs{
        // "a280",
        // "lin318",
        // "d493",
        // "pr1002",
        // "rl1304",
        "d2103",
        "u2319",
        "pcb3038",
        "fnl4461",
        "rl5915",
        "rl5934",
        "pla7397",
    };

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            THEN ("We can generate safe GMI cuts if cutting_loop fails.") {
                OutPrefs prefs;
                Solver solver("problems/" + prob + ".tsp",
                                   //"test_data/tours/" + prob + ".sol",
                                   999, prefs);
                LP::PivType piv = solver.cutting_loop(false, true);

                if (piv == LP::PivType::Frac) {
                    LP::CoreLP &core =
                    const_cast<LP::CoreLP&>(solver.get_core_lp());

                    const vector<double> x = core.lp_vec();
                    const vector<double> tour =
                    solver.tour_basis().best_tour_edges;

                    unique_ptr<Sep::SafeGomory> gmi;

                    Timer t("Safe GMI sep for " + prob);
                    t.start();
                    REQUIRE_NOTHROW(gmi =
                                    util::make_unique<Sep::SafeGomory>(core,
                                                                       tour,
                                                                       x)
                        );

                    bool result;
                    REQUIRE_NOTHROW(result = gmi->find_cuts());
                    CHECK(result);
                    t.stop();
                    t.report(false);
                    cout << "\n";
                }
            }
        }
    }
}

#endif //CMR_HAVE_SAFEGMI
#endif //CMR_DO_TESTS
