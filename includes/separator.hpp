#ifndef CMR_CUTCONTROL_H
#define CMR_CUTCONTROL_H

#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "karp.hpp"

namespace CMR {

/** Namespace for matters related to cuts and separation of cutting planes. */
namespace Sep {

struct CutSelect {
    bool segment = true;
    bool fast2m = true;
    bool exact2m = true;
    bool blkcomb =true ;
    bool chain = false;
    bool gomory = false;
    bool simple_dp = true;
    bool connect_cuts = true;
};

/** Class for separation of cutting planes.
 * This class is instantiated with data about active edges in a relaxation
 * and a current lp solution, then used to search for cuts violated by the
 * solution vector. 
 */
class Separator {
public:

    Separator(Data::GraphGroup &graphdata,
              Data::BestGroup &bestdata,
              Data::SupportGroup &suppdata,
              Data::KarpPartition &kpart);

    Separator(Data::GraphGroup &graphdata,
              Data::BestGroup &bestdata,
              Data::SupportGroup &suppdata,
              Data::KarpPartition &kpart,
              int round_limit);
    
    bool find_cuts(TourGraph &TG);

    bool segment_sep(TourGraph &TG);
    bool fast2m_sep(TourGraph &TG);
    bool blkcomb_sep(TourGraph &TG);

    bool simpleDP_sep();

    bool connect_sep(TourGraph &TG);

    const LPcutList &segment_q() const { return seg_q; }
    const LPcutList &fastblossom_q() const { return fast2m_q; }
    const LPcutList &blockcomb_q() const { return blkcomb_q; }

    const CutQueue<dominoparity> &simpleDP_q() const { return dp_q; }

    const LPcutList &connect_cuts_q() const { return connect_q; }

private:
    const int max_total;
    int running_total;

    CutSelect cut_sel;
    
    Data::GraphGroup &graph_data;
    Data::BestGroup &best_data;
    Data::SupportGroup &supp_data;
    Data::KarpPartition &karp_part;

    std::vector<int> perm_elist;

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CutQueue<dominoparity> dp_q;

    LPcutList connect_q;
};

}
}

#endif
