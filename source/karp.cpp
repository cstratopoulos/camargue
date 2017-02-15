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

KarpPartition::KarpPartition(const Instance &inst) try
{
    int ncount = inst.node_count();
    int seed = inst.seed();
    int norm = inst.ptr()->norm;
    bool dummy_part = false;

    if (ncount < 500) {
        cout << "Problem less than 500 nodes, using trivial Karp partition.\n";
        dummy_part = true;
    }

    if (!dummy_part)
        if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
            dummy_part = true;
            cout << "Incompatible norm, using trivial Karp partition.\n";
        }

    if (dummy_part) {
        part_list.resize(1);
        part_list[0].resize(ncount);
        for (int i = 0; i < ncount; ++i)
            part_list[0][i] = i;
        return;
    }

    CCsubdiv *subdiv_list = nullptr;
    int **part_matrix = nullptr;
    CCrandstate rstate;
    int partcount = 0;

    auto cleanup = util::make_guard([&subdiv_list, &part_matrix, &partcount] {
        CC_IFFREE(subdiv_list, CCsubdiv);
        if (part_matrix) {
            for (int i = 0; i < partcount; ++i)
                CC_IFFREE(part_matrix[i], int);
            CC_FREE(part_matrix, int *);
        }
    });

    CCutil_sprand(seed, &rstate);

    if (CCutil_karp_partition(ncount, inst.ptr(), bucket_size(ncount),
                              &partcount,
                              &subdiv_list, &part_matrix, &rstate))
        throw runtime_error("CCutil_karp_partition failed.");

    part_list.resize(partcount);
    for (int i = 0; i < partcount; ++i) {
        part_list[i].resize(subdiv_list[i].cnt);
        for (int k = 0; k < subdiv_list[i].cnt; ++k)
            part_list[i][k] = part_matrix[i][k];
    }
} catch (const std::exception &e) {
    cerr << e.what() << " in KarpPartition constructor\n";
    throw runtime_error("KarpPartition constructor failed.");
}

int KarpPartition::bucket_size(const int ncount) { return 4 * sqrt(ncount); }

}
}
