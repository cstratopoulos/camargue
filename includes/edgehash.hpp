#ifndef CMR_EDGEHASH_H
#define CMR_EDGEHASH_H

#include "graph.hpp"
#include "util.hpp"

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace CMR {
namespace util {

/** Mapping node pairs to edges. */
class EdgeHash {
public:
    EdgeHash(int size); /**< Construct an EdgeHash for approx size elements. */
    ~EdgeHash();

    EdgeHash(const EdgeHash &eh) = delete;
    EdgeHash &operator=(const EdgeHash &eh) = delete;

    void add(int end1, int end2, int val); /**< Add a pair. */
    void set(int end1, int end2, int val); /**< Set val for existing pair. */
    void erase(int end1, int end2); /**< Delete a pair. */

    std::vector<Edge> get_all(); /**< Get a vector of all the edges. */

    void clear(); /**< Clear all the edges. */

    int get_val(int end1, int end2); /**< Get the val for an edge. */
    

private:
    static constexpr int factor = 1.5;
    CCutil_edgehash eh;
};


}
}

#endif
