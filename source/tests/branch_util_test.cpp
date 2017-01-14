#include "config.hpp"
#include "branch_util.hpp"


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


#ifdef CMR_DO_TESTS

using ObjDig = std::pair<double, int>;

SCENARIO ("Computing branch objective utilities",
          "[ABC][Dive]") {

    vector<ObjDig> test_cases{
        ObjDig(699.0, 3),
        ObjDig(108159.0, 6),
        ObjDig(1234567.0, 7),
    };

    for (ObjDig &p : test_cases) {
        GIVEN ("A double obvjal " + to_string(p.first)) {
            double val = p.first;
            double digits = p.second;

            THEN ("We can compute the number of digits in it") {
                REQUIRE(CMR::ABC::num_digits(val) == digits); 
            }
        }
    }
}

#endif //CMR_DO_TESTS
