#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "edgehash.hpp"

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <utility>

#include <cstdlib>

#include <catch.hpp>

using std::vector;
using std::pair;

using std::string;
using std::to_string;
using std::cout;

SCENARIO ("Instantiating EdgeHash",
          "[util][EdgeHash]") {
    cout << "Sizeof EdgeHash: " << sizeof(CMR::util::EdgeHash) << "\n";
    cout << "Size of CCutil_edgehash: " << sizeof(CCutil_edgehash) << "\n";
    CMR::util::EdgeHash eh(5);

}

#endif //CMR_DO_TESTS
