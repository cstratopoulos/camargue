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
      TG(active_tour_.edges(), core_edges_, active_tour_.tour_perm()),
      perm_elist(s_dat.support_elist)
{
    for (int i = 0; i < perm_elist.size(); ++i)
        perm_elist[i] = active_tour.tour_perm()[perm_elist[i]];
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("MetaCuts constructor failed");
}

bool MetaCuts::find_cuts()
{
    runtime_error err("Problem in MetaCuts::find_cuts");
    double st = util::zeit();

    try {
        if (!attempt_sep())
            return false;
    } CMR_CATCH_PRINT_THROW("testing LP solution", err);

    try {
        if (!price_combs())
            return false;
    } CMR_CATCH_PRINT_THROW("searching for combs within tolerance", err);

    CCtsp_lpgraph lg;
    CC_GCgraph gg;
    int ncount = TG.node_count();
    vector<double> &ecap = supp_data.support_ecap;
    int ecount = ecap.size();

    auto cleanup = util::make_guard([&lg, &gg]
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

    int ic_num = 0;
    for (HGitr it : interest_combs) {
        ++ic_num;
        const HyperGraph &H = *it;
        lpcut_in old;
        lpcut_in *dd = NULL;
        auto lpcut_guard = util::make_guard([&old]
                                            { CCtsp_free_lpcut_in(&old);});

        CCtsp_init_lpcut_in(&old);

        try { old = H.to_lpcut_in(perm); }
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
            lpcut_in *ddnext = dd->next;
            double lp_slack = CCtsp_cutprice(&lg, dd, &ecap[0]);
            double tour_slack = CCtsp_cutprice(TG.pass_ptr(), dd,
                                               TG.tour_array());

            if (lp_slack <= -Epsilon::CutViol &&
                (!filter_primal || tour_slack == 0.0)) {
                found_cuts.emplace_back(dd, lp_slack);
            } else {
                CCtsp_free_lpcut_in(dd);
                CC_IFFREE(dd, CCtsp_lpcut_in);
            }
            dd = ddnext;
        }
    }


    if (found_cuts.empty()) {
        st = util::zeit() - st;
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
        meta_q.push_front(cv.first);

    st = util::zeit() - st;

    if (verbose)
        cout << "\tEnqueued " << meta_q.size() << " " << meta_type
             << " cuts in " << setprecision(2) << st << "s" << setprecision(6)
             << endl;

    return true;
}

bool MetaCuts::price_combs()
{
    runtime_error err("Problem in MetaCuts::price_combs");

    try { price_cliques(); } CMR_CATCH_PRINT_THROW("pricing cliques", err);

    const vector<HyperGraph> &lp_cuts = EC.get_cuts();

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

        if (slack < 10) {
            try { interest_combs.push_back(it); }
            CMR_CATCH_PRINT_THROW("emplacing back comb of interest", err);
        }
    }

    return (!interest_combs.empty());
}

bool MetaCuts::above_threshold(int num_paths)
{
    return num_paths > 0.2 * TG.node_count();
}


/// @returns true iff we should attempt separation. This is always true for
/// Decker cuts. For Handling and Teething it is based on properties of the
/// x-vector.
/// @remark A rewrite of static int no_tighten from concorde/TSP/control.c
bool MetaCuts::attempt_sep()
{
    if (meta_type == Type::Decker)
        return true;

    CC_SRKgraph G;
    auto cleanup = util::make_guard([&G]{ CCcut_SRK_free_graph(&G); });

    int num_paths = 0;
    int ncount = TG.node_count();

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
    if (CCtsp_test_pure_comb(TG.node_count(), &c, &result, (int *) NULL)) {
        cerr << "CCtsp_test_pure_comb failed" << endl;
        throw runtime_error("MetaCuts::pure_comb failed");
    }

    return result;
}

}
}
