#include "karp.hpp"
#include "err_util.hpp"

#include <stdexcept>
#include <iostream>

#include <cmath>

using std::cout;
using std::cerr;

using std::runtime_error;
using std::logic_error;

using std::vector;

namespace CMR {
namespace Data{

KarpPartition::KarpPartition(const int ncount, CCdatagroup *dat,
			     const int seed) try {
  CCsubdiv *subdiv_list = nullptr;
  int **part_matrix = nullptr;
  CCrandstate rstate;
  int partcount = 0;
   
  auto cleanup = CMR::make_guard([&subdiv_list, &part_matrix, &partcount] {
    CC_IFFREE(subdiv_list, CCsubdiv);
    if(part_matrix){
      for(int i = 0; i < partcount; ++i) CC_IFFREE(part_matrix[i], int);
      CC_FREE(part_matrix, int *);
    }
  });

  CCutil_sprand(seed, &rstate);

  if(CCutil_karp_partition(ncount, dat, bucket_size(ncount), &partcount,
			   &subdiv_list, &part_matrix, &rstate))
    throw runtime_error("CCutil_karp_partition failed.");

  part_list.resize(partcount);
  for(int i = 0; i < partcount; ++i){
    part_list[i].resize(subdiv_list[i].cnt);
    for(int k = 0; k < subdiv_list[i].cnt; ++k)
      part_list[i][k] = part_matrix[i][k];
  }
 } catch(const std::exception &e) {
  cerr << e.what() << " in KarpPartition constructor\n";
  throw runtime_error("KarpPartition constructor failed.");
 }

int KarpPartition::bucket_size(const int ncount) { return 2 * sqrt(ncount); }

}
}
