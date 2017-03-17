#include "config.hpp"

#ifdef CMR_DO_TESTS_DISABLED

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

SCENARIO ("Running a Solver with a  DFSbrancher",
          "[ABC][DFSbrancher]") {
    using namespace CMR;
    vector<string> probs{
        "dantzig42",
        "pr76",
        "rat99",
        "a280",
        "lin318",
        "d493",
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

#endif //CMR_DO_TESTS
