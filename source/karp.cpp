#include "karp.hpp"
#include "err_util.hpp"

#include <stdexcept>
#include <string>
#include <iostream>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/util.h>
}

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;

using std::vector;

namespace CMR {
namespace Data{

/**
 * Construct a KarpPartition, with the option to force a dummy partition and/or
 * write the partition to file. A dummy partition consists of the whole vertex
 * set. Dummy partitions are used by default if \p inst has under 300 nodes,
 * or if its norm is incompatible with the Concorde KD tree code.
 * @param inst the Instance to generate the partition from.
 * @param make_dummy if true, force a dummy partition to be used.
 * @param save_part if true, save the partition to a file called
 * `inst.problem_name().kpart`.
 */
KarpPartition::KarpPartition(const Instance &inst, bool make_dummy,
                             bool save_part) try
{
    int ncount = inst.node_count();
    int seed = inst.seed();
    int norm = inst.ptr()->norm;
    bool dummy_part = make_dummy;

    if (make_dummy)
        cout << "Forcing a dummy partition" << endl;

    if (!dummy_part)
        if (ncount < 300) {
            cout << "Problem under 300 nodes, using trivial Karp partition"
                 << endl;
            dummy_part = true;
        }

    if (!dummy_part)
        if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
            dummy_part = true;
            cout << "Incompatible norm, using trivial Karp partition"
                 << endl;
        }

    if (dummy_part) {
        part_list.resize(1);
        part_list[0].resize(ncount);
        for (int i = 0; i < ncount; ++i)
            part_list[0][i] = i;

        if (dummy_part && save_part) {
            cout << "Dummy partition not being written to file." << endl;
        }
        return;
    }

    CCsubdiv *subdiv_list = nullptr;
    int **part_matrix = nullptr;
    CCrandstate rstate;
    int partcount = 0;
    std::string partfile(inst.problem_name() + ".kpart");

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

    if (save_part) {
        if (CCutil_write_subdivision_index(const_cast<char*>(partfile.c_str()),
                                           ncount, partcount, subdiv_list))
            throw runtime_error("CCutil_write_subdivision_index failed.");
        cout << "Wrote partition to file " << partfile <<  ".index " << endl;
    }

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
