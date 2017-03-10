#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "solver.hpp"

#include <catch.hpp>

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

SCENARIO ("Reproducing results with the same random seed",
          "[random_seed]") {
    using namespace CMR;
    vector<string> probs{
        "rat99",
        "lin105",
        "brg180",
        "a280",
        "lin318",
        "d493",
        "p654",
        };

#ifdef CMR_USE_OMP
    throw std::runtime_error("Testing random_seed with OMP enabled");
#endif

    for (string &prob : probs) {
    GIVEN ("Two Solvers for " + prob + " with the same random seed") {
        OutPrefs prefs;
        Solver slv1("problems/" + prob + ".tsp", 99, prefs);
        Solver slv2("problems/" + prob + ".tsp", 99, prefs);

        THEN ("The same cutting loop terminates with the same final results") {
            slv1.cutting_loop(true, true, true);
            slv2.cutting_loop(true, true, true);

            const LP::CoreLP &core1 = slv1.get_core_lp();
            const LP::CoreLP &core2 = slv2.get_core_lp();

            REQUIRE(core1.num_rows() == core2.num_rows());
            REQUIRE(core1.num_cols() == core2.num_cols());
            REQUIRE(core1.get_objval() == core2.get_objval());
        }
    }
    }
}


#endif
