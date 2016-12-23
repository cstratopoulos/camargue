#ifndef CMR_EDGEHASH_H
#define CMR_EDGEHASH_H

#include "util.hpp"

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

namespace CMR {
namespace util {

class EdgeHash {
public:
    EdgeHash(int size);

private:
    static constexpr int factor = 1.5;
    c_struct_ptr<CCutil_edgehash> hash_handle;
};


}
}

#endif
