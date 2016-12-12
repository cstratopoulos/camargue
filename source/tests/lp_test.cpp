#include <catch.hpp>

#include "lp_interface.hpp"
#include "core_lp.hpp"
#include "tests.hpp"
#include "err_util.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::vector;
using std::unique_ptr;
using std::string;
using std::pair;

#ifdef CMR_DO_TESTS

SCENARIO ("Consructing a Core LP",
          "[LP][CoreLP]") {
    vector<string> probs{"dantzig42", "lin318", "pr1002", "pcb3038"};

    for (string& prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("A core LP is constructed") {
                THEN ("Its constructor doesn't throw.") {
                    CMR::Data::Instance inst("problems/" + prob + ".tsp", 99);
                    CMR::Data::GraphGroup g_dat(inst);
                    CMR::Data::BestGroup b_dat(inst, g_dat);
                    REQUIRE_NOTHROW(CMR::LP::CoreLP core(g_dat, b_dat));
                }
            }
        }
    }
}

SCENARIO ("Constructing LP Relaxations",
          "[.LP][.Relaxation][valgrind]") {
    WHEN("A relaxation is constructed"){
        THEN("Its constructor doesn't throw"){
            REQUIRE_NOTHROW(CMR::LP::Relaxation rel);            
        }
    }

    GIVEN ("A constructed Relaxation") {
        CMR::LP::Relaxation rel;
        WHEN ("Another one is constructed") {
            CMR::LP::Relaxation rel2;
            THEN ("It can be move assigned with no leaks or exceptions") {
                REQUIRE_NOTHROW(rel = std::move(rel2));
            }
        }
    }
}

SCENARIO ("Black box testing of failures in constructing LP Relaxations",
          "[.LP][.Relaxaton][!shouldfail][valgrind]") {
    GIVEN ("The raw data structures in an LP Relaxation") {
        WHEN ("A nonzero rval occurs in the constructor") {
            THEN ("No memory is leaked") {
                int rval = 0;
        
                CPXLPptr cplex_lp = (CPXLPptr) NULL;
                CPXENVptr cplex_env = CPXopenCPLEX(&rval);

                rval = 5;

                if (rval) {
                    CPXcloseCPLEX(&cplex_env);
                    REQUIRE_NOTHROW(throw std::runtime_error("rval"));
                }
            }
        }
    }
}

#endif
