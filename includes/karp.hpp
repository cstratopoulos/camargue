/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief KARP PARTITIONS OF TSP INSTANCES
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_KARP_H
#define CMR_KARP_H

#include "datagroups.hpp"

extern "C" {
#include <concorde/INCLUDE/util.h>
}

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
    KarpPartition() = default; /**< Construct an empty partition */

    KarpPartition(const Data::Instance &inst);

  /** Construct a partition from a TSP instance.
   * The data in \p dat is used to construct a Karp partition for a TSP 
   * instance on \p ncount nodes, with random seed \p seed.
   */
  KarpPartition(const int ncount, CCdatagroup *dat, const int seed);

  /** The number of sub-regions into which the data has been partitioned. */
  int num_parts() const { return part_list.size(); }

  /** Returns a vector containing the nodes of the \p i th partition.
   * Valid values are in the range `0...num_parts()`.
   */
  const std::vector<int> &operator[](int i) const { return part_list[i]; }

  /** Maximum partition size as a function of the number of nodes.
   * Reasonable choices are a fraction of \p ncount, or a scalar multiple of
   * `sqrt(ncount)` or something like that.
   * @remark By exhaustively computing simple DP inequalities with no 
   * partitioning, it may be possible to empirically determine the best value
   * for this function.
   */
  static int bucket_size(const int ncount);

private:
  /** A ragged matrix of partitions. */
  std::vector<std::vector<int>> part_list;
};

}
}

#endif
