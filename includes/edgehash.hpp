/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**
 * @file
 * @brief Hash map for graph edges.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_EDGEHASH_H
#define CMR_EDGEHASH_H

#include "graph.hpp"
#include "util.hpp"

#include <iostream>
#include <memory>
#include <utility>
#include <vector>

namespace CMR {
namespace util {

/// Hash map for node pairs representing edges.
/// @remark A wrapper to Concorde's CCutil_edgehash.
class EdgeHash {
public:
    EdgeHash(int size); //!< An EdgeHash for approximately \p size elements.
    ~EdgeHash();

    EdgeHash(const EdgeHash &eh) = delete; //!< No copy constructor.
    EdgeHash &operator=(const EdgeHash &eh) = delete; //!< No copy assign.

    void add(int end1, int end2, int val); //!< Add a pair.
    void set(int end1, int end2, int val); //!< Set val for existing pair.
    void erase(int end1, int end2); //!< Delete a pair.

    std::vector<Graph::Edge> get_all(); //!< Get a vector of all the edges.
    void clear(); //!< Clear all the edges from the hash.
    int get_val(int end1, int end2); //!< Get the val for an edge.


private:
    struct eh_impl; //!< The hash table implementation.
    std::unique_ptr<eh_impl> eh_pimpl; //!< Pimpl idiom firewall.
};


}
}

#endif
