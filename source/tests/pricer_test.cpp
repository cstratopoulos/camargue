#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "solver.hpp"
#include "util.hpp"
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



SCENARIO ("Running the Solver cutting_loop on augmentable or optimal instances",
          "[Price][Pricer][add_edges]") {
    vector<string> probs{
        "pr76",
        "lin105",
        "pcb442",
        "p654",
    };

    for (string &prob : probs) {
        GIVEN (prob) {
            THEN ("We can run cutting loop and price on opt/aug tours") {
                CMR::OutPrefs prefs;
                CMR::Solver solver("problems/" + prob + ".tsp",
                                   //prob + ".sol",
                                   0, prefs);

                REQUIRE_NOTHROW(solver.cutting_loop());
            }
        }
    }
}

#endif //CMR_DO_TESTS
