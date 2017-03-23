#include "pool_sep.hpp"
#include "err_util.hpp"
#include "timer.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

using std::cout;
using std::cerr;
using std::endl;
using std::setprecision;

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::vector;

namespace CMR {

namespace Eps = Epsilon;

namespace Sep {

using CutType = HyperGraph::Type;
using DoublePair = std::pair<double, double>;


PoolCuts::PoolCuts(ExternalCuts &EC_,
                   const std::vector<Graph::Edge> &core_edges_,
                   const LP::ActiveTour &active_tour_,
                   Data::SupportGroup &s_dat) try
    : EC(EC_), core_edges(core_edges_), active_tour(active_tour_),
      tour_edges(active_tour_.edges()),
      supp_data(s_dat), hg_q(50),
      lp_slacks(vector<double>(EC_.get_cutpool().size(), 0.0)),
      tour_slacks(vector<double>(EC_.get_cutpool().size(), 0.0))
{
    clique_vals.reserve(EC.get_pool_cbank().size());

    vector<int> keep_inds;
    int ncount = supp_data.supp_graph.node_count;

    keep_inds.reserve(ncount);
    for (int i = 0; i < tour_edges.size(); ++i)
        if (tour_edges[i] > 1.0 - Epsilon::Zero)
            keep_inds.push_back(i);

    tour_adj = Graph::AdjList(ncount, core_edges, tour_edges, keep_inds);
    if (tour_adj.edge_count != ncount) {
        cerr << "Tour adj edge count " << tour_adj.edge_count << endl;
        throw runtime_error("Not all edges made it into tour adj");
    }

} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("PoolCuts constructor failed.");
}

bool PoolCuts::tighten_pool()
{
    double st = util::zeit();
    using lpcut_in = CCtsp_lpcut_in;
    runtime_error err("Problem in PoolCuts::tighten_pool");

    try {
        if (!price_cuts(true))
            return false;
    } CMR_CATCH_PRINT_THROW("pricing cuts", err);

    TourGraph TG;
    try {
        TG = TourGraph(tour_edges, core_edges, active_tour.tour_perm());
    } CMR_CATCH_PRINT_THROW("constructing TourGraph", err);

    vector<int> lp_perm_elist;
    const vector<int> &perm = active_tour.tour_perm();

    try {
        lp_perm_elist = supp_data.support_elist;
        for (int i = 0; i < lp_perm_elist.size(); ++i)
            lp_perm_elist[i] = perm[lp_perm_elist[i]];
    } CMR_CATCH_PRINT_THROW("building lp perm elist", err);

    CCtsp_lpgraph lg;
    int ncount = supp_data.supp_graph.node_count;
    vector<double> &ecap = supp_data.support_ecap;
    int ecount = ecap.size();

    auto cleanup = util::make_guard([&lg] { CCtsp_free_lpgraph (&lg); });

    CCtsp_init_lpgraph_struct(&lg);

    if (CCtsp_build_lpgraph(&lg, ncount, ecount, &lp_perm_elist[0],
                            (int *) NULL)) {
        cerr << "CCtsp_build_lpgraph failed" << endl;
        throw err;
    }

    if (CCtsp_build_lpadj(&lg, 0, ecount)) {
        cerr << "CCtsp_build_lpadj failed" << endl;
        throw err;
    }

    CCtsp_tighten_info tstats; //unused
    CCtsp_init_tighten_info(&tstats);

    vector<HyperGraph> &cutpool = EC.cut_pool;
    int prev_size = cutpool.size();
    using PoolItr = vector<HyperGraph>::iterator;

    // for cuts of interest, store their CC pointer, lp slack, and iterator
    // to cut they were derived from.
    using CutViolPitr = std::tuple<lpcut_in *, double, PoolItr>;
    vector<CutViolPitr> found_cuts;

    for (int i = 0; i < cutpool.size(); ++i) {
        if (!slack_of_interest(lp_slacks[i], tour_slacks[i], true))
            continue;
        double lp_slack = lp_slacks[i];
        double tour_slack = tour_slacks[i];
        HyperGraph &H = cutpool[i];

        lpcut_in old;
        lpcut_in *tight = NULL;
        bool want_new = false;
        auto lpcut_guard = util::make_guard([&old, &want_new, &tight]
                                            {
                                                CCtsp_free_lpcut_in(&old);
                                                if (!want_new) {
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

        double lp_act = 0.0;

        /* Split cases based on filter_primal and lp_slack/tour_slack */
        /* If filter_primal...
               - If the cut is already tight at the tour, try to tighten it
                 wrt the lp, then check it is still tight at tour.
               - Else, try to tighten it wrt the tour, then check if it is
               violated by the LP. */
        /* Else....
              - Just try and tighten it wrt the LP. */

        if (filter_primal) {
            if (tour_slack == 0.0) { //already tight, tigthen wrt lp
                double lp_improve = 0.0;
                if (CCtsp_tighten_lpcut_in(&lg, &old, &ecap[0], tight,
                                           &tstats, &lp_improve)) {
                    cerr << "CCtsp_tighten_lpcut_in failed" << endl;
                    throw err;
                }
                lp_act = lp_slack - lp_improve;

                // check still tight at tour
                if (lp_act <= -Epsilon::CutViol) {
                    double tour_slack = CCtsp_cutprice(TG.pass_ptr(), tight,
                                                       TG.tour_array());
                    want_new = (tour_slack == 0.0);
                }
            } else { //tighten wrt tour
                double tour_improve = 0.0;
                if (CCtsp_tighten_lpcut_in(TG.pass_ptr(), &old,
                                           TG.tour_array(), tight, &tstats,
                                           &tour_improve)) {
                    cerr << "CCtsp_tighten_lpcut_in failed" << endl;
                    throw err;
                }
                double tour_act = tour_slack - tour_improve;

                // want it if violated
                if (tour_act == 0.0) {
                    cout << "MADE A SLACK PRIMAL CUT TIGHT!!!!!" << endl;
                    lp_act = CCtsp_cutprice(&lg, tight, &ecap[0]);
                    want_new = (lp_act <= -Epsilon::CutViol);
                    if (want_new)
                        cout << "AND LP VIOLATED!!!!!" << endl;
                }
            }
        } else { // don't care abt primal, tighten wrt lp
            double lp_improve = 0.0;
            if (CCtsp_tighten_lpcut_in(&lg, &old, &ecap[0], tight,
                                       &tstats, &lp_improve)) {
                cerr << "CCtsp_tighten_lpcut_in failed" << endl;
                throw err;
            }
            lp_act = lp_slack - lp_improve;
            want_new = (lp_act <= -Epsilon::CutViol);
        }

        if (want_new)
            try {
                found_cuts.emplace_back(tight, lp_act, cutpool.begin() + i);
            } CMR_CATCH_PRINT_THROW("emplacing successful tighten", err);
    }

    if (found_cuts.empty()) {
        st = util::zeit() - st;
        if (verbose)
            cout << "\tFound 0 tighten_pool cuts in "
                 << setprecision(2) << st << "s" << setprecision(6)
                 << "(" << prev_size << " cuts in pool)" << endl;
        return false;
    }

    if (found_cuts.size() > 250) {
        std::sort(found_cuts.begin(), found_cuts.end(),
                  [](const CutViolPitr &A, const CutViolPitr &B)
                  { return std::get<1>(A) < std::get<1>(B); });

        while(found_cuts.size() > 250) {
            lpcut_in *delp = std::get<0>(found_cuts.back());
            CCtsp_free_lpcut_in(delp);
            CC_IFFREE(delp, CCtsp_lpcut_in);
            found_cuts.pop_back();
        }
    }

    for (const CutViolPitr &CVP : found_cuts) {
        HyperGraph &H = *std::get<2>(CVP);
        H = HyperGraph();

        lpcut_in *found = std::get<0>(CVP);
        tight_q.push_front(found);
    }

    cutpool.erase(std::remove_if(cutpool.begin(), cutpool.end(),
                                 [](const HyperGraph &H)
                                 { return H.cut_type() == CutType::Non; }),
                  cutpool.end());

    int after_size = cutpool.size();

    st = util::zeit() - st;

    if (verbose)
        cout << "\tEnqueued " << tight_q.size() << " tighten pool cuts in "
             << setprecision(2) << st << "s" << setprecision(6)
             << " (Cutpool size " << prev_size << " -> "
             << after_size << ")" << endl;

    return true;
}

bool PoolCuts::find_cuts()
{
    runtime_error err("Problem in PoolCuts::find_cuts");

    if (verbose)
        cout << "Cutpool sep, filter_primal " << filter_primal << ", "
             << EC.cut_pool.size() << " cuts in pool." << endl;

    try {
        if (!price_cuts(false))
            return false;
    } CMR_CATCH_PRINT_THROW("pricing cuts", err);

    vector<HyperGraph> &cutpool = EC.cut_pool;

    using PoolItr = vector<HyperGraph>::iterator;
    using ItrIdx = std::pair<PoolItr, int>;
    vector<ItrIdx> primal_cuts;

    try {
        for (int i = 0; i < cutpool.size(); ++i)
            if (slack_of_interest(lp_slacks[i], tour_slacks[i], false))
                primal_cuts.push_back(ItrIdx(cutpool.begin() + i, i));
    } CMR_CATCH_PRINT_THROW("creating iterator vector", err);

    int qcap = hg_q.q_capacity();

    if (primal_cuts.size() > qcap) {
        std::sort(primal_cuts.begin(), primal_cuts.end(),
                  [this](const ItrIdx &I, const ItrIdx &J)
                  { return lp_slacks[I.second] < lp_slacks[J.second]; });
        primal_cuts.resize(qcap);
    }

    for (ItrIdx &I : primal_cuts) {
        hg_q.emplace_back(std::move(*(I.first)));
        *(I.first) = HyperGraph();
    }

    if (verbose)
        cout << "\t" << hg_q.size() << " cuts removed from pool for addition."
             << endl;

    cutpool.erase(std::remove_if(cutpool.begin(), cutpool.end(),
                                 [](const HyperGraph &H)
                                 { return H.cut_type() == CutType::Non; }),
                  cutpool.end());

    return true;
}

bool PoolCuts::price_cuts(bool tighten)
{
    Timer t("Price cuts");
    t.start();

    Timer pc("Price cliques", &t);
    pc.start();
    price_cliques();
    pc.stop();

    Timer hgt("Pricing HG cuts", &t);
    Timer dpt("Pricing dominos", &t);

    const vector<HyperGraph> &pool = EC.get_cutpool();
    bool result = false;
    int numfound = 0;

    for (int i = 0; i < pool.size(); ++i) {
        const HyperGraph &H = pool[i];
        CutType h_type = H.cut_type();
        if (h_type == CutType::Non || h_type == CutType::Branch)
            throw logic_error("Pool has an invalid cut");

        double lp_slack = 0.0;
        double tour_slack = 0.0;
        double rhs = H.get_rhs();

        if (h_type == CutType::Comb) {
            hgt.resume();
            lp_slack -= rhs;
            tour_slack -= rhs;
            for (const Clique::Ptr &clq_ref : H.get_cliques()) {
                DoublePair p = clique_vals[*clq_ref];
                lp_slack += p.first;
                tour_slack += p.second;
            }
            hgt.stop();
        } else if (h_type == CutType::Domino) { // we just expand the whole row
            if (!tighten) {
                dpt.resume();
                LP::SparseRow R;

                try {
                    H.get_coeffs(core_edges, R.rmatind, R.rmatval);
                } catch (const exception &e) {
                    cerr << e.what() << "\n";
                    throw runtime_error("Couldn't get DP row");
                }

                lp_slack = rhs - get_activity(supp_data.lp_vec, R);
                tour_slack = rhs - get_activity(tour_edges, R);

                dpt.stop();
            } else {
                lp_slack = 1000;
                tour_slack = 1000;
            }
        }

        lp_slacks[i] = lp_slack;
        tour_slacks[i] = tour_slack;
        if (slack_of_interest(lp_slack, tour_slack, tighten)) {
            result = true;
            ++numfound;
        }
    }

    t.stop();

    if (verbose > 1) {
        cout << "\tPriced cuts, found " << numfound << " of interest" << endl;
            pc.report(false);
            hgt.report(false);
            dpt.report(false);
            t.report(false);
    }
    
    return result;
}


void PoolCuts::price_cliques()
{
    vector<Graph::Node> &lp_nodelist = supp_data.supp_graph.nodelist;
    vector<Graph::Node> &tour_nodelist = tour_adj.nodelist;

    for (Graph::Node &n : lp_nodelist)
        n.mark = 0;

    const CliqueBank &pool_cliques = EC.get_pool_cbank();
    const vector<int> &def_tour = pool_cliques.ref_tour();

    int marker = 0;

    for (CliqueBank::ConstItr it = pool_cliques.begin();
         it != pool_cliques.end(); ++it) {
        const Clique &clq = it->first;
        double lp_val = 0.0;
        double tour_val = 0.0;

        ++marker;

        for (const Segment &seg : clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                lp_nodelist[node].mark = marker;
                tour_nodelist[node].mark = marker;
            }

        for (const Segment &seg : clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];

                for (Graph::AdjObj &a : lp_nodelist[node].neighbors)
                    if (lp_nodelist[a.other_end].mark != marker)
                        lp_val += a.val;

                for (Graph::AdjObj &a : tour_nodelist[node].neighbors)
                    if (tour_nodelist[a.other_end].mark != marker)
                        tour_val += a.val;
            }
        clique_vals[clq] = DoublePair(lp_val, tour_val);
    }
}

inline bool PoolCuts::slack_of_interest(double lp_slack, double tour_slack,
                                        bool tighten)
{
    if (tighten)
        return (lp_slack < 0.5 && (!filter_primal || tour_slack < 0.5));
    else
        return (lp_slack <= -Eps::CutViol &&
                (!filter_primal || tour_slack == 0.0));
}

}
}
