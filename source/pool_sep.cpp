#include "pool_sep.hpp"
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


PoolCuts::PoolCuts(ExternalCuts &EC_,
                   const std::vector<Graph::Edge> &core_edges_,
                   Data::SupportGroup &s_dat) try
    : EC(EC_), core_edges(core_edges_), supp_data(s_dat),
      hg_q(50), cut_slacks(vector<double>(EC_.get_cutpool().size(), 0.0))
{
    clique_vals.reserve(EC.get_pool_cbank().size());
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("PoolCuts constructor failed.");
}

bool PoolCuts::find_cuts()
{
    runtime_error err("Problem in PoolCuts::find_cuts");
    bool result = false;

    try { result = price_cuts(); } CMR_CATCH_PRINT_THROW("pricing cuts", err);
    if (!result)
        return result;

    using CVpair = std::pair<HyperGraph, double>;
    vector<CVpair> viol_cuts;
    vector<HyperGraph> &cut_pool = EC.cut_pool;

    try {
        for (int i = 0; i < cut_pool.size(); ++i)
            if (cut_slacks[i] <= -Eps::Cut) {
                viol_cuts.emplace_back(CVpair(std::move(cut_pool[i]),
                                              cut_slacks[i]));
                cut_pool[i] = HyperGraph();
            }
    } CMR_CATCH_PRINT_THROW("emplacing violated cut pairs", err);

    cut_pool.erase(std::remove_if(cut_pool.begin(), cut_pool.end(),
                                  [](const HyperGraph &H)
                                  { return H.cut_type() == CutType::Non; }),
                   cut_pool.end());

    int end_range = 0;
    int q_cap = hg_q.q_capacity();

    if (viol_cuts.size() > q_cap) {
        end_range = q_cap;
        std::sort(viol_cuts.begin(), viol_cuts.end(),
                  [](const CVpair &A, const CVpair &B)
                  { return A.second > B.second; }); //most violated cuts at back
    } else
        end_range = viol_cuts.size();

    cout << "Viol cuts size " << viol_cuts.size() << endl;
    cout << "End range " << end_range << endl;

    for (auto r_it = viol_cuts.rbegin(); r_it != r_it + end_range; ++r_it) {
        hg_q.emplace_back(std::move(r_it->first));
        viol_cuts.pop_back();
    }

    try {
        for (CVpair &CV : viol_cuts)
            cut_pool.emplace_back(std::move(CV.first));
    } CMR_CATCH_PRINT_THROW("putting cuts back in pool", err);

    return result;
}

bool PoolCuts::price_cuts()
{
    bool report = true;
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

    for (int i = 0; i < pool.size(); ++i) {
        const HyperGraph &H = pool[i];
        CutType h_type = H.cut_type();
        if (h_type == CutType::Non || h_type == CutType::Branch)
            throw logic_error("Pool has an invalid cut");

        double slack = 0.0;
        double rhs = H.get_rhs();

        if (h_type == CutType::Comb) {
            hgt.resume();
            slack -= rhs;
            for (const Clique::Ptr &clq_ref : H.get_cliques())
                slack += clique_vals[*clq_ref];
            hgt.stop();
        } else if (h_type == CutType::Domino) { // we just expand the whole row
            dpt.resume();
            LP::SparseRow R;

            try {
                H.get_coeffs(core_edges, R.rmatind, R.rmatval);
            } catch (const exception &e) {
                cerr << e.what() << "\n";
                throw runtime_error("Couldn't get DP row");
            }

            double activity = get_activity(supp_data.lp_vec, R);
            slack = rhs - activity;
            dpt.stop();
        }

        cut_slacks[i] = slack;
        if (slack <= -Eps::Cut)
            result = true;
    }

    t.stop();

    if (report) {
        pc.report(false);
        hgt.report(false);
        dpt.report(false);
        t.report(false);
    }


    return result;
}


void PoolCuts::price_cliques()
{
    vector<Graph::Node> &nodelist = supp_data.supp_graph.nodelist;

    for (Graph::Node &n : nodelist)
        n.mark = 0;

    const CliqueBank &pool_cliques = EC.get_pool_cbank();
    const vector<int> &def_tour = pool_cliques.ref_tour();
    int marker = 0;

    for (CliqueBank::ConstItr it = pool_cliques.begin();
         it != pool_cliques.end(); ++it) {
        const Clique &clq = it->first;
        double val = 0.0;
        ++marker;
        for (const Segment &seg : clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                nodelist[node].mark = marker;
            }

        for (const Segment &seg : clq.seg_list())
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                for (Graph::AdjObj &a : nodelist[node].neighbors)
                    if (nodelist[a.other_end].mark != marker)
                        val += a.val;
            }
        clique_vals[clq] = val;
    }
}

}
}
