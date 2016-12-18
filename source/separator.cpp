#include "separator.hpp"
#include "cc_lpcuts.hpp"
#include "simpleDP.hpp"
#include "err_util.hpp"

#include <vector>
#include <stdexcept>
#include <iostream>

using std::vector;

using std::runtime_error;
using std::logic_error;

using std::cout;


namespace CMR {
namespace Sep {

bool Separator::find_cuts(CMR::TourGraph &TG)
{
    runtime_error err("Problem in Separator::find_cuts.");
    
    vector<int> perm_elist;
    vector<double> &ecap = supp_data.support_ecap;
    vector<int> &perm = best_data.perm;

    try {
        perm_elist = supp_data.support_elist;
    } CMR_CATCH_PRINT_THROW("allocating permuted elist", err);

    for(int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = perm[perm_elist[i]];

    int running_total = 0;
    bool found_segments = false;
    bool found_comb = false;

    try {
        SegmentCuts segments(perm_elist, ecap, TG, seg_q);

        if (segments.find_cuts()) {
            found_segments = true;
//            cout << "\t"  << seg_q.size() << " segment cuts\n";
            running_total += seg_q.size();
            if (running_total > 8)
                return true;
        }

        FastBlossoms fast2m(perm_elist, ecap, TG, fast2m_q);

        if (fast2m.find_cuts()) {
            found_comb = true;
//            cout << "\t"  << fast2m_q.size() << " fast blossoms\n";
            running_total += fast2m_q.size();
            if (running_total > 8)
                return true;
        }

        BlockCombs blkcomb(perm_elist, ecap, TG, blkcomb_q);

        if (blkcomb.find_cuts()) {
            found_comb = true;
//            cout << "\t"  << blkcomb_q.size() << " block combs\n";
            running_total += blkcomb_q.size();
            if (running_total > 8)
                return true;
        }

        if (supp_data.connected && !found_segments && !found_comb) {
            if (supp_data.in_subtour_poly()) {
                dp_q.clear();
                SimpleDP dominos(graph_data, karp_part,
                                 best_data, supp_data, dp_q);
                if (dominos.find_cuts()) {
                    cout << "\t" << dp_q.size() << " dominos\n";
                    running_total += dp_q.size();
                }
            }
        }
        
    } CMR_CATCH_PRINT_THROW("calling primal sep", err);

    if (running_total > 0)
        return true;

    if (!supp_data.connected && supp_data.integral) {
        try {
            ConnectCuts subtours(perm_elist, ecap, TG, connect_q);
            
            if (subtours.find_cuts()) {
//                cout << "\t" << connect_q.size() << " connect cuts\n";
                running_total += connect_q.size();
            }
        } CMR_CATCH_PRINT_THROW("calling subtour sep", err);
    }

    return (running_total > 0);
}

}
}
