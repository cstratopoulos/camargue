#ifndef CMR_CUTCONTROL_H
#define CMR_CUTCONTROL_H

#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "karp.hpp"

namespace CMR {

/** Namespace for matters related to cuts and separation of cutting planes. */
namespace Sep {

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

    LPcutList seg_q;
    LPcutList fast2m_q;
    LPcutList blkcomb_q;

    CutQueue<dominoparity> dp_q;

    LPcutList connect_q;

    friend class Solver;

private:
    const int max_total;
    int running_total;
    
    Data::GraphGroup &graph_data;
    Data::BestGroup &best_data;
    Data::SupportGroup &supp_data;
    Data::KarpPartition &karp_part;

    std::vector<int> perm_elist;
};

}
}

#endif
