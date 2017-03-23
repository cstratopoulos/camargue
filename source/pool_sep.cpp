#include "pool_sep.hpp"
#include "cc_lpcuts.hpp"
#include "err_util.hpp"
#include "timer.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>

using std::cout;
using std::cerr;
using std::endl;

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

// bool PoolCuts::tighten_pool()
// {
//     runtime_error err("Problem in PoolCuts::tighten_pool");

//     if (verbose)
//         cout << "Cutpool tighten sep, filter_primal " << filter_primal
//              << ", " << EC.cut_pool.size() << " cuts in pool." << endl;

//     TourGraph TG;
//     try {
//         TG = TourGraph(tour_edges, )
//     }

//     vector<int> tour_perm_elist;
//     vector<int> lp_perm_elist;

//     try {
//         if (!price_cuts(true))
//             return true;
//     } CMR_CATCH_PRINT_THROW("pricing cuts", err);
// }

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
            if (lp_slacks[i] <= -Eps::CutViol &&
                (tour_slacks[i] == 0 || filter_primal == false))
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

                double lp_activity = get_activity(supp_data.lp_vec, R);
                lp_slack = rhs - lp_activity;

                if (filter_primal)
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

    if (verbose) {
        cout << "\tPriced cuts, found" << numfound << endl;
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
