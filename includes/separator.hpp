#ifndef CMR_CUTCONTROL_H
#define CMR_CUTCONTROL_H

#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "karp.hpp"

#include <memory>

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
    Separator(Data::GraphGroup &graphdata,
              Data::BestGroup &bestdata,
              Data::SupportGroup &suppdata,
              Data::KarpPartition &kpart,
              Graph::TourGraph &_TG);

    /// Construct a Separator, but with limit of number of cuts returned.
    Separator(Data::GraphGroup &graphdata,
              Data::BestGroup &bestdata,
              Data::SupportGroup &suppdata,
              Data::KarpPartition &kpart,
              Graph::TourGraph &_TG,
              int round_limit);

    bool segment_sep();
    bool fast2m_sep();
    bool blkcomb_sep();
    bool exact2m_sep();

    bool simpleDP_sep();

    bool connect_sep();

    const LPcutList &segment_q() const { return seg_q; }
    const LPcutList &fastblossom_q() const { return fast2m_q; }
    const LPcutList &blockcomb_q() const { return blkcomb_q; }

    const CutQueue<ex_blossom> &exblossom_q() const { return ex2m_q; }
    const CutQueue<dominoparity> &simpleDP_q() const { return dp_q; }

    const LPcutList &connect_cuts_q() const { return connect_q; }

private:
    int max_total;
    int running_total;

    Data::GraphGroup &graph_data;
    Data::BestGroup &best_data;
    Data::SupportGroup &supp_data;
    Data::KarpPartition &karp_part;

    Graph::TourGraph &TG;

    std::vector<int> perm_elist;

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CutQueue<ex_blossom> ex2m_q;
    CutQueue<dominoparity> dp_q;

    LPcutList connect_q;
};

/// This is a shorthand for resetting a unique_ptr to a Separator, like
/// a surrogate move assign since Separator has mostly reference members.
void ptr_reset(std::unique_ptr<Separator> &sep, Data::GraphGroup &g_dat,
               Data::BestGroup &b_dat, Data::SupportGroup &s_dat,
               Data::KarpPartition &kpart,
               Graph::TourGraph &TG);


}
}

#endif
