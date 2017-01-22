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

/// Mapping node pairs to edges via wrapper to Concorde data structure.
class EdgeHash {
public:
    /// Construct an Edgehash for approximately \p size elements.
    EdgeHash(int size);
    ~EdgeHash(); //!< Destruct and free associated memory.

    EdgeHash(const EdgeHash &eh) = delete; //!< No copy constructor.
    EdgeHash &operator=(const EdgeHash &eh) = delete; //!< No copy assign.

    void add(int end1, int end2, int val); //!< Add a pair.
    void set(int end1, int end2, int val); //!< Set val for existing pair.
    void erase(int end1, int end2); //!< Delete a pair.

    std::vector<Graph::Edge> get_all(); //!< Get a vector of all the edges.
    void clear(); //!< Clear all the edges from the hash. 
    int get_val(int end1, int end2); //!< Get the val for an edge.
    

private:
    /// A scaling factor used to set the capacity of the underlying array.
    static constexpr int factor = 1.5;
    CCutil_edgehash eh; //!< The underlying Concorde data structure.
};


}
}

#endif
