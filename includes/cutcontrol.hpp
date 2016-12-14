#ifndef CMR_CUTCONTROL_H
#define CMR_CUTCONTROL_H

#include "cuts.hpp"
#include "core_lp.hpp"
#include "datagroups.hpp"
#include "karp.hpp"

namespace CMR {
namespace Sep {

class CutControl {
public:
    CutControl(CMR::Data::GraphGroup &graphdata,
               CMR::Data::BestGroup &bestdata,
               CMR::Data::SupportGroup &suppdata,
               CMR::Data::KarpPartition &kpart) :
        graph_data(graphdata), best_data(bestdata), supp_data(suppdata),
        karp_part(kpart){}

    bool find_cuts(CMR::TourGraph &TG);

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CMR::CutQueue<CMR::dominoparity> dp_q;

    LPcutList connect_q;

    const CMR::Data::SupportGroup &support_data() const { return supp_data; }

private:
    CMR::Data::GraphGroup &graph_data;
    CMR::Data::BestGroup &best_data;
    CMR::Data::SupportGroup &supp_data;
    CMR::Data::KarpPartition &karp_part;
};

}
}

#endif
