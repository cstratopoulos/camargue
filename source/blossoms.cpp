#include "blossoms.hpp"
#include "config.hpp"
#include "err_util.hpp"
#include "util.hpp"

extern "C" {
  #include <concorde/INCLUDE/cut.h>
}

#include <iostream>
#include <vector>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;

using std::vector;

using std::runtime_error;
using std::exception;

namespace CMR {
namespace Eps = Epsilon;
namespace Sep {

ExBlossoms::ExBlossoms(const vector<Graph::Edge> &_edges,
                       const vector<double> &x,
                       Data::BestGroup &b_dat, Data::SupportGroup &s_dat,
                       CutQueue<ex_blossom> &_blossom_q) :
    edges(_edges), lp_vec(x),
    best_data(b_dat), supp_data(s_dat), blossom_q(_blossom_q) {}

bool ExBlossoms::find_cuts() {
    runtime_error err("Problem in ExBlossoms::find_cuts");

    vector<int> &sup_inds = supp_data.support_indices;
    vector<double> &sup_ecap = supp_data.support_ecap;
    vector<int> &sup_elist = supp_data.support_elist;
    const vector<int> &tour_edges = best_data.best_tour_edges;
    
    vector<double> cut_ecap;

    try { cut_ecap = sup_ecap; }
    CMR_CATCH_PRINT_THROW("copying cut ecap", err);

    for (auto i = 0; i < sup_inds.size(); ++i) {
        int index = sup_inds[i];
        if (tour_edges[index] == 1)
            cut_ecap[i] = 1 - sup_ecap[i];
    }

    vector<ex_blossom> intermediate_cuts;
    int ncount = best_data.best_tour_nodes.size();

    for (auto i = 0; i < sup_inds.size(); ++i) {
        int cut_ind = sup_inds[i];
        int tour_entry = tour_edges[cut_ind];
        EndPts e(sup_elist[2 * i], sup_elist[(2 * i) + 1]);

        int cut_count = 0;
        int *cut_nodes = (int *) NULL;
        double orig_weight = 1.0;
        double changed_weight = 1.0;
        double cut_val = 1.0;

        if (tour_entry == 0) {
            orig_weight = sup_ecap[i];
            changed_weight = 1 - sup_ecap[i];
        } else if (tour_entry == 1) {
            orig_weight = 1 - sup_ecap[i];
            changed_weight = sup_ecap[i];
        }

        cut_ecap[i] = changed_weight;

        //reverts cut_ecap when it goes out of scope.
        auto ecap_guard = util::make_guard([&cut_ecap, i, orig_weight]
                                           { cut_ecap[i] = orig_weight; });

        if (CCcut_mincut_st(ncount, sup_inds.size(),
                            &sup_elist[0], &cut_ecap[0],
                            e.end[0], e.end[1],
                            &cut_val, &cut_nodes, &cut_count)) {
            cerr << "CCcut_mincut_st failed.\n";
            throw err;
        }

        //frees cut nodes when it goes out of scope
        util::c_array_ptr<int> cnodes_ptr(cut_nodes);

        if (cut_val >= 1.0 - Eps::Cut ||
            cut_count < 3 || cut_count > ncount - 3)
            continue;

        vector<int> handle;
        
        try {
            for (int j = 0; j < cut_count; ++j)
                handle.push_back(cut_nodes[j]);
        } CMR_CATCH_PRINT_THROW("copying handle", err);

        try {
            intermediate_cuts.emplace_back(handle, cut_ind, cut_val);
        } CMR_CATCH_PRINT_THROW("emplacing intermediate cut", err);
    }

    if (intermediate_cuts.empty())
        return false;

    try {
        intermediate_cuts.erase(std::remove_if(intermediate_cuts.begin(),
                                               intermediate_cuts.end(),
                                               [this, &tour_edges, ncount]
                                               (const ex_blossom &B)
                                               {
                                                   return
                                                   bad_blossom(B, tour_edges,
                                                               lp_vec, edges,
                                                               ncount);
                                               }),
                                intermediate_cuts.end());
    } CMR_CATCH_PRINT_THROW("filtering bad blossoms", err);

    if (intermediate_cuts.empty()) {
        cout << "\n|||Cuts were found but had to be cleaned up|||\n\n";
        return false;
    }
    
    std::sort(intermediate_cuts.begin(), intermediate_cuts.end(),
              [](const ex_blossom &B, const ex_blossom &C)
              { return B.cut_val < C.cut_val; });

    for (auto it = intermediate_cuts.rbegin(); it != intermediate_cuts.rend();
         ++it) {
        blossom_q.push_front(*it);
    }

    return true;
}

}
}
