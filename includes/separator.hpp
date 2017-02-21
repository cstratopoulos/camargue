#ifndef CMR_CUTCONTROL_H
#define CMR_CUTCONTROL_H

#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "karp.hpp"
#include "graph.hpp"

#include <memory>
#include <vector>

namespace CMR {

/// Matters related to cuts and separation of cutting planes.
namespace Sep {

/** Class for separation of cutting planes.
 * This class is instantiated with data about active edges in a relaxation
 * and a current lp solution, then used to search for cuts violated by the
 * solution vector.
 */
class Separator {
public:

    /// Construct a Separator using problem data.
    Separator(const std::vector<Graph::Edge> &core_edges_,
              Data::BestGroup &bestdata,
              Data::SupportGroup &suppdata,
              Data::KarpPartition &kpart);

    bool segment_sep();
    bool fast2m_sep();
    bool blkcomb_sep();
    bool exact2m_sep();

    bool simpleDP_sep();

    bool connect_sep();

    LPcutList &segment_q()  { return seg_q; }
    LPcutList &fastblossom_q()  { return fast2m_q; }
    LPcutList &blockcomb_q()  { return blkcomb_q; }

    CutQueue<ex_blossom> &exblossom_q()  { return ex2m_q; }
    CutQueue<dominoparity> &simpleDP_q()  { return dp_q; }

    LPcutList &connect_cuts_q()  { return connect_q; }

private:
    const std::vector<Graph::Edge> &core_edges;
    Data::BestGroup &best_data;
    Data::SupportGroup &supp_data;
    Data::KarpPartition &karp_part;

    TourGraph TG;

    std::vector<int> perm_elist;

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CutQueue<ex_blossom> ex2m_q;
    CutQueue<dominoparity> dp_q;

    LPcutList connect_q;
};

}
}

#endif
