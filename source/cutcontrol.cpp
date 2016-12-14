#include "cutcontrol.hpp"
#include "cc_lpcuts.hpp"
#include "simpleDP.hpp"
#include "err_util.hpp"

#include <vector>
#include <stdexcept>

using std::vector;

using std::runtime_error;
using std::logic_error;


namespace CMR {
namespace Sep {

bool CutControl::find_cuts(CMR::TourGraph &TG)
{
    runtime_error err("Problem in CutControl::find_cuts.");
    
    vector<int> perm_elist;
    vector<double> &ecap = supp_data.support_ecap;

    try {
        perm_elist = supp_data.support_elist;
    } CMR_CATCH_PRINT_THROW("allocating permuted elist", err);

    for (int &i : perm_elist)
        i = best_data.perm[i];

    int running_total = 0;
    bool found_segments = false;

    try {
        SegmentCuts segments(perm_elist, ecap, TG, seg_q);

        if (segments.find_cuts()) {
            found_segments = true;
            running_total += seg_q.size();
            if (running_total > 15)
                return true;
        }

        FastBlossoms fast2m(perm_elist, ecap, TG, fast2m_q);

        if (fast2m.find_cuts()) {
            running_total += fast2m_q.size();
            if (running_total > 15)
                return true;
        }

        BlockCombs blkcomb(perm_elist, ecap, TG, blkcomb_q);

        if (blkcomb.find_cuts()) {
            running_total += blkcomb_q.size();
            if (running_total > 15)
                return true;
        }

        if (supp_data.connected && !found_segments) {
            if (supp_data.in_subtour_poly()) {
                dp_q.clear();
                SimpleDP dominos(graph_data, karp_part,
                                 best_data, supp_data, dp_q);
                if (dominos.find_cuts())
                    running_total += dp_q.size();
            }
        }
        
    } CMR_CATCH_PRINT_THROW("calling primal sep", err);

    if (running_total > 0)
        return true;

    if (!supp_data.connected && supp_data.integral) {
        try {
            ConnectCuts subtours(perm_elist, ecap, TG, connect_q);
            
            if (subtours.find_cuts())
                running_total += connect_q.size();
        } CMR_CATCH_PRINT_THROW("calling subtour sep", err);
    }

    return (running_total > 0);
}

}
}
