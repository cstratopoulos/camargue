/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Some straightforward separation routines.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_SEPARATOR_H
#define CMR_SEPARATOR_H

#include "active_tour.hpp"
#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "karp.hpp"
#include "graph.hpp"
#include "hypergraph.hpp"

#include <memory>
#include <vector>

namespace CMR {

/// Matters related to cuts and separation of cutting planes.
namespace Sep {

/// Management of basic separation routines.
class Separator {
public:
    /// Construct a Separator using problem data.
    Separator(const std::vector<Graph::Edge> &core_edges_,
              const LP::ActiveTour &active_tour_,
              Data::SupportGroup &suppdata,
              Data::KarpPartition &kpart, int seed);

    bool segment_sep();
    bool fast2m_sep();
    bool blkcomb_sep();
    bool exact2m_sep();

    bool simpleDP_sep();

    bool connect_sep();
    bool exsub_sep();

    bool pool_sep(ExternalCuts &EC);
    bool tighten_pool(ExternalCuts &EC);
    bool consec1_sep(ExternalCuts &EC);

    bool local_sep(int chunk_sz, bool sphere);

    LPcutList &segment_q()  { return seg_q; }
    LPcutList &fastblossom_q()  { return fast2m_q; }
    LPcutList &blockcomb_q()  { return blkcomb_q; }

    CutQueue<ex_blossom> &exblossom_q()  { return ex2m_q; }
    CutQueue<dominoparity> &simpleDP_q()  { return dp_q; }

    LPcutList &connect_cuts_q()  { return connect_q; }
    LPcutList &exact_sub_q() { return exsub_q; }

    LPcutList &cutpool_q() { return pool_q; }
    LPcutList &consec1_q() { return con1_q; }

    LPcutList &local_cuts_q() { return local_q; }

    /// The desired filter_primal value for ConcordeSeparator inheritors.
    bool filter_primal = true;

    bool verbose = false;

private:
    void set_TG(); //!< Construct the TourGraph TG.

    const std::vector<Graph::Edge> &core_edges;
    const LP::ActiveTour &active_tour;
    Data::SupportGroup &supp_data;
    Data::KarpPartition &karp_part;

    TourGraph TG;

    std::vector<int> perm_elist;

    const int random_seed;

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CutQueue<ex_blossom> ex2m_q;
    CutQueue<dominoparity> dp_q;

    LPcutList connect_q;
    LPcutList exsub_q;

    LPcutList pool_q;
    LPcutList con1_q;

    LPcutList local_q;
};

}
}

#endif
