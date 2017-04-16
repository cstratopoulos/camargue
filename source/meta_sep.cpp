#include "meta_sep.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::setprecision;

using std::exception;
using std::runtime_error;
using std::logic_error;

using std::vector;

namespace CMR {
namespace Sep {

using lpcut_in = CCtsp_lpcut_in;

std::ostream &operator<<(std::ostream &os, MetaCuts::Type t)
{
    using Mtype = MetaCuts::Type;
    if (t == Mtype::Decker)
        os << "Double Decker";
    else if (t == Mtype::Handling)
        os << "Handling";
    else if (t == Mtype::Teething)
        os << "Teething";
    else
        throw std::logic_error("Unimplemented MetaCuts::Type <<");

    return os;
}

MetaCuts::MetaCuts(const ExternalCuts &EC_,
                   const vector<Graph::Edge> &core_edges_,
                   const LP::ActiveTour &active_tour_,
                   Data::SupportGroup &s_dat) try
    : EC(EC_), core_edges(core_edges_), active_tour(active_tour_),
      supp_data(s_dat),
      TG(),
      perm_elist(s_dat.support_elist)
{
    if (filter_primal)
        TG = TourGraph(active_tour.edges(), core_edges,
                       active_tour.tour_perm());
    for (int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = active_tour.tour_perm()[perm_elist[i]];
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("MetaCuts constructor failed");
}

bool MetaCuts::tighten_cuts()
{
    double st = util::zeit();

    runtime_error err("Problem in MetaCuts::tighten_cuts");

    try {
        if (!price_combs(true))
            return false;
    } CMR_CATCH_PRINT_THROW("searching for combs within tolerance", err);

    CCtsp_lpgraph lg;
    int ncount = supp_data.supp_graph.node_count;
    vector<double> &ecap = supp_data.support_ecap;
    int ecount = ecap.size();

    auto cleanup = util::make_guard([&lg] { CCtsp_free_lpgraph (&lg); });

    CCtsp_init_lpgraph_struct(&lg);

    if (CCtsp_build_lpgraph(&lg, ncount, ecount, &perm_elist[0],
                            (int *) NULL)) {
        cerr << "CCtsp_build_lpgraph failed" << endl;
        throw err;
    }

    if (CCtsp_build_lpadj(&lg, 0, ecount)) {
        cerr << "CCtsp_build_lpadj failed" << endl;
        throw err;
    }

    using CutViol = std::pair<lpcut_in *, double>;
    vector<CutViol> found_cuts;
    const std::vector<int> &perm = active_tour.tour_perm();

    CCtsp_tighten_info tstats; // unused
    CCtsp_init_tighten_info(&tstats);

    int ic_num = -1;
    for (HGitr it : interest_combs) {
        ++ic_num;
        const HyperGraph &H = *it;
        lpcut_in old;
        lpcut_in *tight = NULL;
        bool want_cut = false;
        auto cutsguard = util::make_guard([&old, &tight, &want_cut]
                                          {
                                              CCtsp_free_lpcut_in(&old);
                                              if (!want_cut) {
                                                  CCtsp_free_lpcut_in(tight);
                                                  CC_IFFREE(tight,
                                                            CCtsp_lpcut_in);
                                              }
                                          });

        tight = CC_SAFE_MALLOC(1, CCtsp_lpcut_in);
        if (tight == NULL) {
            cerr << "Couldn't allocate tight" << endl;
            throw err;
        }

        CCtsp_init_lpcut_in(&old);
        CCtsp_init_lpcut_in(tight);

        try { old = H.to_lpcut_in(perm, true); }
        CMR_CATCH_PRINT_THROW("getting lpcut_in/skelly from HyperGraph", err);

        double improve = 0.0;
        if (CCtsp_tighten_lpcut_in(&lg, &old, &ecap[0], tight, &tstats,
                                   &improve)) {
            cerr << "CCtsp_tighten_lpcut_in failed" << endl;
            throw err;
        }

        double lp_slack = 0.0;
        double tour_slack = 0.0;

        if (improve != 0) {
            lp_slack = CCtsp_cutprice(&lg, tight, &ecap[0]);
            if (lp_slack <= -Epsilon::CutViol) {
                if (filter_primal)
                    tour_slack = CCtsp_cutprice(TG.pass_ptr(), tight,
                                                TG.tour_array());
                if (tour_slack == 0.0) {
                    want_cut = true;
                    try { found_cuts.emplace_back(tight, lp_slack); }
                    CMR_CATCH_PRINT_THROW("emplacing cut", err);
                }
            }
        }
    }

    if (found_cuts.empty()) {
        st = util::zeit() - st;
        if (verbose)
            cout << "\tFound 0 Tighten LP cuts in "
                 << setprecision(2) << st << "s" << setprecision(6) << endl;
        return false;
    }

    if (found_cuts.size() > 250) {
        std::sort(found_cuts.begin(), found_cuts.end(),
                  [](const CutViol &a, const CutViol &b)
                  { return a.second < b.second; });

        while (found_cuts.size() > 250) {
            CCtsp_free_lpcut_in(found_cuts.back().first);
            CC_IFFREE(found_cuts.back().first, CCtsp_lpcut_in);
            found_cuts.pop_back();
        }
    }

    for (const CutViol &cv : found_cuts)
        meta_q.push(cv.first);

    st = util::zeit() - st;

    if (verbose)
        cout << "\tEnqueued " << meta_q.size() << " Tighten LP cuts in "
             << setprecision(2) << st << "s" << setprecision(6) << endl;

    return true;
}

bool MetaCuts::find_cuts(Type meta_type)
{
    runtime_error err("Problem in MetaCuts::find_cuts");
    double st = util::zeit();

    try { //always attempt for Decker cuts
        if (meta_type != Type::Decker && !attempt_sep())
                return false;
    } CMR_CATCH_PRINT_THROW("testing LP solution", err);

    try {
        if (!price_combs(false))
            return false;
    } CMR_CATCH_PRINT_THROW("searching for combs within tolerance", err);

    CCtsp_lpgraph lg;
    CC_GCgraph gg;
    int ncount = supp_data.supp_graph.node_count;
    vector<double> &ecap = supp_data.support_ecap;
    int ecount = ecap.size();

    auto graph_cleanup = util::make_guard([&lg, &gg]
                                          {
                                              CCtsp_free_lpgraph (&lg);
                                              CCcombs_GC_free_graph (&gg);
                                          });

    CCtsp_init_lpgraph_struct(&lg);
    CCcombs_GC_init_graph(&gg);

    if (CCtsp_build_lpgraph(&lg, ncount, ecount, &perm_elist[0],
                            (int *) NULL)) {
        cerr << "CCtsp_build_lpgraph failed" << endl;
        throw err;
    }

    if (CCtsp_build_lpadj(&lg, 0, ecount)) {
        cerr << "CCtsp_build_lpadj failed" << endl;
        throw err;
    }

    if (meta_type == Type::Decker || meta_type == Type::Handling)
        if (CCcombs_GC_build_graph(&gg, ncount, ecount, &perm_elist[0],
                                   &ecap[0])) {
            cerr << "CCcombs_GC_build_graph failed" << endl;
            throw err;
        }

    using CutViol = std::pair<lpcut_in *, double>;
    vector<CutViol> found_cuts;
    const std::vector<int> &perm = active_tour.tour_perm();

    int ic_num = -1;
    for (HGitr it : interest_combs) {
        ++ic_num;
        const HyperGraph &H = *it;
        lpcut_in old;
        lpcut_in *dd = NULL;
        auto lpcut_guard = util::make_guard([&old]
                                            { CCtsp_free_lpcut_in(&old);});

        CCtsp_init_lpcut_in(&old);

        try { old = H.to_lpcut_in(perm, false); }
        CMR_CATCH_PRINT_THROW("getting lpcut_in from HyperGraph", err);

        try {
            if (!pure_comb(old))
                continue;
        } CMR_CATCH_PRINT_THROW("testing old cut", err);

        switch(meta_type) {
        case Type::Decker:
            if (CCtsp_comb_to_double_decker(&lg, &gg, &ecap[0], &old, &dd)) {
                cerr << "CCtsp_comb_to_double_decker failed" << endl;
                throw err;
            }
            break;
        case Type::Handling:
            if (CCtsp_comb_handling(&lg, &gg, &ecap[0], &old, &dd)) {
                cerr << "CCtsp_comb_handling failed" << endl;
                throw err;
            }
            break;
        case Type::Teething:
            if (CCtsp_teething(&lg, &ecap[0], &old, &dd)) {
                cerr << "CCtsp_teething failed" << endl;
                throw err;
            }
            break;
        default:
            throw logic_error("Unimplemented caller in MetaCuts::find_cuts");
        }

        while (dd) {
            lpcut_in *old_dd = dd;
            lpcut_in *ddnext = dd->next;
            double lp_slack = CCtsp_cutprice(&lg, dd, &ecap[0]);
            double tour_slack = 0.0;
            bool want_cut = false;

            auto ddg = util::make_guard([&old_dd, &want_cut]
                                        {
                                            if (!want_cut) {
                                                CCtsp_free_lpcut_in(old_dd);
                                                CC_IFFREE(old_dd,
                                                          CCtsp_lpcut_in);
                                            }
                                        });

            if (lp_slack <= -Epsilon::CutViol) {
                if (filter_primal)
                    tour_slack = CCtsp_cutprice(TG.pass_ptr(), dd,
                                                TG.tour_array());
                if (tour_slack == 0.0) {
                    want_cut = true;
                    try { found_cuts.emplace_back(dd, lp_slack); }
                    CMR_CATCH_PRINT_THROW("emplacing cut", err);
                }
            }
            dd = ddnext;
        }
    }


    if (found_cuts.empty()) {
        st = util::zeit() - st;
        if (verbose)
            cout << "\tFound 0 " << meta_type << " cuts in "
                 << setprecision(2) << st << "s" << setprecision(6) << endl;
        return false;
    }

    if (found_cuts.size() > 250) {
        std::sort(found_cuts.begin(), found_cuts.end(),
                  [](const CutViol &a, const CutViol &b)
                  { return a.second < b.second; });

        while (found_cuts.size() > 250) {
            CCtsp_free_lpcut_in(found_cuts.back().first);
            CC_IFFREE(found_cuts.back().first, CCtsp_lpcut_in);
            found_cuts.pop_back();
        }
    }

    for (const CutViol &cv : found_cuts)
        meta_q.push(cv.first);

    st = util::zeit() - st;

    if (verbose)
        cout << "\tEnqueued " << meta_q.size() << " " << meta_type
             << " cuts in " << setprecision(2) << st << "s" << setprecision(6)
             << endl;

    return true;
}

bool MetaCuts::price_combs(bool tighten)
{
    runtime_error err("Problem in MetaCuts::price_combs");

    try { price_cliques(); } CMR_CATCH_PRINT_THROW("pricing cliques", err);

    const vector<HyperGraph> &lp_cuts = EC.get_cuts();

    double tolerance = tighten ? 0.5 : 10;

    for (HGitr it = lp_cuts.begin(); it != lp_cuts.end(); ++it) {
        const HyperGraph &H = *it;
        int num_cliques = H.get_cliques().size();
        double rhs = H.get_rhs();

        if (H.cut_type() != HyperGraph::Type::Comb ||
            num_cliques < 4 || (num_cliques % 2) != 0 || H.get_sense() != 'G')
            continue;

        double slack = -rhs;
        for (const Clique::Ptr &clq_ref : H.get_cliques())
            slack += clique_vals[*clq_ref];

        if (slack < tolerance) {
            try { interest_combs.push_back(it); }
            CMR_CATCH_PRINT_THROW("emplacing back comb of interest", err);
        }
    }

    return (!interest_combs.empty());
}

bool MetaCuts::above_threshold(int num_paths)
{
    return num_paths > 0.2 * supp_data.supp_graph.node_count;
}


/// @returns true iff we should attempt separation. This is always true for
/// Decker cuts. For Handling and Teething it is based on properties of the
/// x-vector.
/// @remark A rewrite of static int no_tighten from concorde/TSP/control.c
bool MetaCuts::attempt_sep()
{
    CC_SRKgraph G;
    auto cleanup = util::make_guard([&G]{ CCcut_SRK_free_graph(&G); });

    int num_paths = 0;
    int ncount = supp_data.supp_graph.node_count;

    CCcut_SRK_init_graph(&G);
    if (CCcut_SRK_buildgraph(&G, ncount, supp_data.support_ecap.size(),
                             &supp_data.support_elist[0],
                             &supp_data.support_ecap[0]))
        throw runtime_error("CCcut_SRK_buildgraph failed");

    CCcut_SRK_increment_marker(&G);

    if (CCcut_SRK_defluff(&G))
        throw runtime_error("CCcut_SRK_defluff failed");

    CCcut_SRK_identify_paths_to_edges(&G, &num_paths, 0);

    return above_threshold(num_paths);
}


void MetaCuts::price_cliques()
{
    try { clique_vals.reserve(EC.get_cbank().size()); }
    catch (const exception &e) {
        cerr << e.what() << " reserving clique_vals" << endl;
        throw runtime_error("MetaCuts::price_cliques failed");
    }

    vector<Graph::Node> &lp_nodelist = supp_data.supp_graph.nodelist;

    for (Graph::Node &n : lp_nodelist)
        n.mark = 0;

    const CliqueBank &lp_cliques = EC.get_cbank();
    const vector<int> &def_tour = lp_cliques.ref_tour();

    int marker = 0;

    for (CliqueBank::ConstItr it = lp_cliques.begin();
         it != lp_cliques.end(); ++it) {
        const Clique &clq = it->first;
        double lp_val = 0.0;

        ++marker;

        for (const Segment &seg : clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                lp_nodelist[node].mark = marker;
            }

        for (const Segment &seg : clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];

                for (Graph::AdjObj &a : lp_nodelist[node].neighbors)
                    if (lp_nodelist[a.other_end].mark != marker)
                        lp_val += a.val;
            }
        clique_vals[clq] = lp_val;
    }
}

bool MetaCuts::pure_comb(lpcut_in &c)
{
    int result = 0;
    int ncount = supp_data.supp_graph.node_count;
    if (CCtsp_test_pure_comb(ncount, &c, &result, (int *) NULL)) {
        cerr << "CCtsp_test_pure_comb failed" << endl;
        throw runtime_error("MetaCuts::pure_comb failed");
    }

    return result;
}

}
}
