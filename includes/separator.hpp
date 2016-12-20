#ifndef CMR_CUTCONTROL_H
#define CMR_CUTCONTROL_H

#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "karp.hpp"

namespace CMR {

namespace Sep {

/** Class for separation of cutting planes.
 * This class is instantiated with data about active edges in a relaxation
 * and a current lp solution, then used to search for cuts violated by the
 * solution vector. 
 */
class Separator {
public:

    Separator(CMR::Data::GraphGroup &graphdata,
              CMR::Data::BestGroup &bestdata,
              CMR::Data::SupportGroup &suppdata,
              CMR::Data::KarpPartition &kpart) :
        graph_data(graphdata), best_data(bestdata), supp_data(suppdata),
        karp_part(kpart){}
    
    bool find_cuts(CMR::TourGraph &TG);

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CMR::Sep::CutQueue<CMR::Sep::dominoparity> dp_q;

    LPcutList connect_q;

    const CMR::Data::SupportGroup &support_data() const { return supp_data; }

    friend class Solver;

private:
    CMR::Data::GraphGroup &graph_data;
    CMR::Data::BestGroup &best_data;
    CMR::Data::SupportGroup &supp_data;
    CMR::Data::KarpPartition &karp_part;
};

}
}

#endif
