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

PoolCuts::PoolCuts(ExternalCuts &EC_,
                   const std::vector<Graph::Edge> &core_edges_,
                   Data::SupportGroup &s_dat) try
    : EC(EC_), core_edges(core_edges_), supp_data(s_dat),
      hg_q(50), cut_slacks(vector<double>(EC_.get_cutpool().size(), 0.0))
{
    clique_vals.reserve(EC.get_cbank().size());
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

    return result;
}

bool PoolCuts::price_cuts()
{
    Timer t("Price cuts");
    t.start();

    Timer pc("Price cliques", &t);
    pc.start();
    price_cliques();
    pc.stop();

    Timer hgt("Pricing HG cuts", &t);
    Timer dpt("Pricing dominos", &t);

    using CutType = HyperGraph::Type;
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
        if (slack <= Eps::Cut)
            result = true;
    }

    t.stop();

    pc.report(false);
    hgt.report(false);
    dpt.report(false);
    t.report(false);


    return result;
}


void PoolCuts::price_cliques()
{
    vector<Graph::Node> &nodelist = supp_data.supp_graph.nodelist;

    for (Graph::Node &n : nodelist)
        n.mark = 0;

    const CliqueBank &cbank = EC.get_cbank();
    const vector<int> &def_tour = cbank.ref_tour();
    int marker = 0;

    for (CliqueBank::ConstItr it = cbank.begin(); it != cbank.end(); ++it) {
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
