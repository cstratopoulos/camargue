#include <catch.hpp>

#include "lp_interface.hpp"
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

SCENARIO ("Constructing LP Relaxations",
          "[LP][Relaxation]") {
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
