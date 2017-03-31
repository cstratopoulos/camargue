/**
 * @file
 * @brief Pricer implementations for exact lower bounds/edge elimination.
 */

#include "pricer.hpp"

#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <utility>

using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::exception;
using std::runtime_error;
using std::logic_error;

namespace CMR {

using LP::PivType;
using CutType = Sep::HyperGraph::Type;
using f64 = util::Fixed64;

namespace Eps = Epsilon;

namespace Price {

f64 Pricer::exact_lb(bool full)
{
    vector<PrEdge<f64>> priced_edges;
    return exact_lb(full, priced_edges);
}

f64 Pricer::exact_lb(bool full,
                     vector<PrEdge<f64>> &priced_edges)
{
    for (const Sep::HyperGraph &H : ext_cuts.get_cuts())
        if (H.cut_type() == CutType::Non)
            throw logic_error("Pricer::exact_lb was called w non cut present");

    runtime_error err("Problem in Pricer::exact_lb");

    try {
        ex_duals = util::make_unique<LP::DualGroup<f64>>(true, core_lp,
                                                         ext_cuts);
    } CMR_CATCH_PRINT_THROW("constructing exact DualGroup", err);

    priced_edges.clear();

    vector<double> d_pi;
    vector<char> senses;
    vector<double> rhs_vec;

    int numrows = core_lp.num_rows();
    int numcols = core_lp.num_cols();
    int ncount = inst.node_count();

    vector<f64> ex_pi;

    try {
        core_lp.get_pi(d_pi, 0, numrows - 1);
        senses = core_lp.senses(0, numrows - 1);
        core_lp.get_rhs(rhs_vec, 0, numrows -1);

        ex_pi = vector<f64>(d_pi.begin(), d_pi.end());
    } CMR_CATCH_PRINT_THROW("grabbing cut info", err);

    for (int i = ncount; i < numrows; ++i)
        if (senses[i] == 'G') {
            if (ex_pi[i] < 0.0)
                ex_pi[i] = 0.0;
            else if (senses[i] == 'L')
                if (ex_pi[i] < 0.0)
                    ex_pi[i] = 0.0;
        }

    f64 pi_sum{0.0};

    for (int i = 0; i < numrows; ++i)
        util::add_mult(pi_sum, ex_pi[i], rhs_vec[i]);


    vector<PrEdge<f64>> target_edges;

    if (full) {
        for (int i = 0; i < ncount; ++i)
            for (int j = 0; j < ncount; ++j)
                target_edges.emplace_back(i, j);
    } else {
        for (const Graph::Edge &e : core_graph.get_edges())
            target_edges.emplace_back(e.end[0], e.end[1]);
    }

    price_edges(target_edges, ex_duals, true);

    f64 rc_sum{0.0};

    for (const PrEdge<f64> &e : target_edges)
        if (e.redcost < 0.0)
            rc_sum -= e.redcost;

    f64 bound = pi_sum - rc_sum;

    priced_edges = std::move(target_edges);

    return bound;
}

void Pricer::elim_edges(bool make_opt)
{
    runtime_error err("Problem in Pricer::elim_edges.");

    if (make_opt) {
        double ot = util::zeit();
        cout << "Optimizing before elimination...";
        try { core_lp.primal_opt(); } CMR_CATCH_PRINT_THROW("optimizing", err);
        cout << "obj val " << core_lp.get_objval() << " in "
             << (util::zeit() - ot) << "s" << endl;
    }

    f64 tourlen{core_lp.global_ub()};
    f64 lower_bd{0.0};
    vector<PrEdge<f64>> graph_edges;

    try {
        lower_bd = exact_lb(false, graph_edges);
    } CMR_CATCH_PRINT_THROW("getting exact lb", err);

    f64 gap{tourlen - lower_bd};
    f64 cutoff{gap - 1};

    if (cutoff < 0) {
        if (verbose)
            cout << "Negative cutoff, do not elim." << endl;
        return;
    }

    if (verbose)
        cout << "\tElimination cutoff " << cutoff << endl;

    try {
        util::ptr_reset(ex_duals, true, core_lp, core_lp.external_cuts());
    } CMR_CATCH_PRINT_THROW("getting exact duals", err);

    vector<int> col_delset;
    int ecount = core_graph.edge_count();

    try {
        col_delset.resize(ecount, 0);
    } CMR_CATCH_PRINT_THROW("reserving/prepping vectors", err);

    const vector<int> &tour_colstat = core_lp.get_active_tour().base().colstat;
    const vector<double> &tour_edges = core_lp.get_active_tour().edges();

    int elimct = 0;

    for (int i = 0; i < ecount; ++i) {
        const PrEdge<f64> &e = graph_edges[i];
        if (tour_edges[i] != 0.0 || tour_colstat[i] != 0)
            continue;

        if (e.redcost > cutoff) {
            col_delset[i] = 1;
            ++elimct;
        }
    }

    cout << "\t\t" << elimct << " edges can be eliminated\n" << endl;
    if (elimct == 0)
        return;

    try {
        core_lp.remove_edges(std::move(col_delset));
    } CMR_CATCH_PRINT_THROW("adjusting CoreLP/CoreGraph stuff", err);
}

bool Pricer::scan_edges(vector<PrEdge<f64>> &gen_edges, int &loop1,
                        int &loop2)
{
    if (!ex_duals)
        throw logic_error("Tried to scan edges without exact duals.");

    const vector<f64> &node_pi_est = ex_duals->node_pi_est;
    const Graph::AdjList &alist = core_graph.get_adj();

    int ncount = inst.node_count();
    int i = loop1;
    int j = loop2;
    int first = 1;

    gen_edges.clear();

    if (i >= ncount)
        return true;

    for (; i < ncount; ++i) {
        for (const Graph::AdjObj &a : alist.nodelist[i].neighbors) {
            int j = a.other_end;
            if (j > i) {
                f64 len = a.val;
                f64 rc = len - node_pi_est[i] - node_pi_est[j];
                if (rc < 0.0)
                    gen_edges.emplace_back(i, j, rc);
                if (gen_edges.size() == f64Batch) {
                    loop1 = i;
                    return false;
                }
            }
        }
    }

    // for(; i < ncount; ++i) {
    //     int stop = ncount;
    //     if (first == 0)
    //         j = i + 1;
    //     first = 0;
    //     for(; j < stop; ++j) {
    //         int end = j;
    //         f64 rc = inst.edgelen(i, j) - node_pi_est[i] - node_pi_est[j];
    //         if (rc < 0.0) {
    //             gen_edges.emplace_back(i, end, rc);
    //             if (gen_edges.size() == f64Batch) {
    //                 loop1 = i;
    //                 loop2 = j + 1;
    //                 return false;
    //             }
    //         }
    //     }
    // }

    loop1 = ncount;
    loop2 = ncount;
    return true;
}

}
}
