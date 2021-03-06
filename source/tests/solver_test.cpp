#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "solver.hpp"
#include "util.hpp"
#include "timer.hpp"


#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>

#include <cstdlib>

#include <catch.hpp>

using std::min;
using std::max;
using std::array;
using std::vector;

using std::string;
using std::to_string;
using std::cout;

SCENARIO ("Running the Solver cutting_loop with no pricing.",
          "[Solver][cutting_loop]") {
    vector<string> probs{
        "dantzig42",
        "pr76",
        "lin318",
        "d493",
        "p654",
        "pr1002",
        "d2103",
        "pr2392"
    };

    for (string &prob : probs) {
        GIVEN (prob) {
            THEN ("We can instantiate a Solver and run cutting_loop") {
                CMR::OutPrefs prefs;
                CMR::Solver solver("problems/" + prob + ".tsp",
                                   //"test_data/tours/" + prob + ".sol",
                                   99, prefs);

                REQUIRE_NOTHROW(solver.cutting_loop(false, true, true));
            }
        }
    }
}

SCENARIO ("Running random 1k solver cutting_loop",
          "[Solver][cutting_loop][random]") {
    for (int i : {1, 2, 3, 4, 5}) {
        GIVEN ("A 1k random instance" + to_string(i)) {
            THEN ("We can instantiate a Solver and run cutting_loop") {
                CMR::OutPrefs prefs;
                CMR::Solver solver(0, 1000, 1000000, prefs);

                REQUIRE_NOTHROW(solver.cutting_loop(false, true, true));
            }
        }
    }
}

#endif //CMR_DO_TESTS
