#ifndef PSEP_KARP_H
#define PSEP_KARP_H

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <vector>

namespace PSEP {
namespace Data {

class KarpPartition {
public:
  KarpPartition() = default;
  KarpPartition(const int ncount, CCdatagroup *dat, const int seed);

  int num_parts() const { return part_list.size(); }

  const std::vector<int> &operator[](int i) const { return part_list[i]; }

private:
  std::vector<std::vector<int>> part_list;
  static int bucket_size(const int ncount);
};

}
}

#endif
