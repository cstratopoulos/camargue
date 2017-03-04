/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Karp partitions of TSP instances.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_KARP_H
#define CMR_KARP_H

#include "datagroups.hpp"

#include <vector>

namespace CMR {
namespace Data {

/** Class for computing and storing Karp partitions.
 * Computes a geometric partition of a TSP instance, a la Probabilistic
 * Analysis of Partitioning Algorithms for the Traveling-Salesman Problem in
 * The Plane (Karp, 1977). Given geometric TSP data, it partitions the set
 * of nodes into rectangles containing at most some fixed number of nodes.
 * @remark This class is used exclusively in the separation
 * of simple domino parity inequalities, but it shall be owned by a
 * GraphGroup because it only needs to be computed once for a given instance
 * and it is invariant under changes to the active edge set/support graph.
 * @remark This is a highly limited interface to the Concorde function
 * CCutil_karp_partition and the Concorde structure CCsubdiv.
 */
class KarpPartition {
public:
    KarpPartition() = default; //!< Default construct an empty partition.

    KarpPartition(const Data::Instance &inst, bool make_dummy,
                  bool save_part);

    /// Construct a partition from an Instance \p inst.
    KarpPartition(const Data::Instance &inst)
        : KarpPartition(inst, false, false) {}


    /// The number of sub-regions into which the data has been partitioned.
    int num_parts() const { return part_list.size(); }

    const std::vector<int> &operator[](int i) const
        { return part_list[i]; } //!< The nodes of the \p i th partition.

    /// Maximum partition size as function of number of nodes \p ncount.
    static int bucket_size(const int ncount);

private:
    /// A ragged matrix storing the partitions.
    std::vector<std::vector<int>> part_list;
};

}
}

#endif
