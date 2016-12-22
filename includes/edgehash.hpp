#ifndef CMR_EDGEHASH_H
#define CMR_EDGEHASH_H

#include "util.hpp"

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <iostream>
#include <stdexcept>

namespace CMR {
namespace util {

template <int size>
class EdgeHash {
public:
    EdgeHash();
    
    
    ~EdgeHash();

private:
    c_struct_ptr<CCutil_edgehash> hash_handle;
};

template<int size>
EdgeHash<size>::EdgeHash() try
{
    CCutil_edgehash *eh = new CCutil_edgehash;

    hash_handle(eh, CCutil_edgehash_free);
    
} catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    throw std::runtime_error("EdgeHash constructor failed.");
}

}
}

#endif
