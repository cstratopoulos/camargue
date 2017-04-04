/**
 * @file
 * @brief Private Solver methods with loop calls.
 */

#include "config.hpp"
#include "solver.hpp"
#include "separator.hpp"

#if CMR_HAVE_SAFEGMI
#include "safeGMI.hpp"
#endif

#include "err_util.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <functional>
#include <vector>

#include <cmath>
#include <cstdio>

using std::abs;
using std::ceil;

using std::function;
using std::vector;
using std::unique_ptr;

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

using CutType = Sep::HyperGraph::Type;
using PivType = LP::PivType;
namespace Eps = Epsilon;

/**
 * This function template is used to add one class of cuts at a time, attempt
 * to measure progress obtained, and return information on whether a further
 * separation routine should be called.
 * @tparam Qtype the queue representation of the cuts found by \p sepcall.
 * @param[in] sepcall a function which returns true if cuts are found, in which
 * case they are stored in \p sep_q.
 * @param[out] piv if cuts are found, this is the pivot type that occurred
 * after they were added.
 * @param[in,out] piv_stats tracker for the pivot objective values seen this
 * round.
 */
template<typename Qtype>
bool Solver::call_separator(const function<bool()> &sepcall, Qtype &sep_q,
                            PivType &piv, PivStats &piv_stats)
{
    bool result = sepcall();
    if (result) {
        core_lp.pivot_back(true);
        core_lp.add_cuts(sep_q);
        piv = core_lp.primal_pivot();

        double prev_val = piv_stats.prev_val;
        double new_val = core_lp.get_objval();
        double tourlen = core_lp.active_tourlen();

        piv_stats.update(new_val, tourlen);
        piv_stats.found_cuts = true;

        if (output_prefs.prog_bar)
            place_pivot(0.99 * tourlen, tourlen, new_val);

        if (output_prefs.verbose) {
            printf("\t^^Pivot %.2f -> %.2f, ratio %.2f\n",
                   prev_val, new_val, piv_stats.delta_ratio);
            piv_stats.report_extrema();
        }
    }

    return result;
}

double Solver::PH_delta(double new_val, double prev_val, double tourlen)
{
    return std::abs((new_val - prev_val) / (tourlen - prev_val));
}

void Solver::PivStats::update(double new_val, double tourlen)
{
    delta_ratio = Solver::PH_delta(new_val, prev_val, tourlen);
    first_last_ratio = Solver::PH_delta(new_val, initial_piv, tourlen);
    prev_val = new_val;

    if (delta_ratio > max_ratio)
        max_ratio = delta_ratio;

    if (prev_val < lowest_piv)
        lowest_piv = prev_val;
}

void Solver::PivStats::report_extrema()
{
    printf("\tLowest piv %.2f\tMax ratio %.2f\tFirst-last ratio %.2f\n",
           lowest_piv, max_ratio, first_last_ratio);
    cout << std::flush;
}

/**
 * @param piv the PivType of the last call to primal_pivot.
 * @param delta_metric some measurement of the change obtained from the last
 * round of cut generation.
 */
bool Solver::restart_loop(LP::PivType piv, double delta_metric)
{
    return (piv == PivType::Subtour || !core_lp.supp_data.connected ||
            delta_metric >= Eps::PHratio);
}

inline bool Solver::return_pivot(LP::PivType piv)
{
    return (piv == PivType::Tour || piv == PivType::FathomedTour);
}

inline void Solver::place_pivot(double low_lim, double best_tourlen,
                                double piv_val)
{
    int target_entry = 1;
    char piv_char = '>';
    if (piv_val <= low_lim)
        piv_char = '!';
    else if (piv_val == best_tourlen) {
        piv_char = '*';
        target_entry = 78;
    } else {
        target_entry = ceil(78 * ((piv_val - low_lim) /
                                  (best_tourlen - low_lim)));
    }

    p_bar[target_entry] = piv_char;

    int i = 0;
    cout << p_bar[i];
    for (i = 1; i < target_entry; ++i)
        cout << '#';
    while (i < 80)
        cout << p_bar[i++];

    cout << "\r" << std::flush;
    p_bar[target_entry] = ' ';
}

/**@name Convenience macros for separation routines in Solver::cut_and_piv.
 * These macros are meant to be used in the body of the main `while` loop
 * in Solver::cut_and_piv to invoke a separation routine using
 * Solver::call_separator, and then restart or advance the loop as needed.
 * This saves the use of repetitive boilerplate which can make it hard to
 * navigate the loop at a glance, given the `return` and `continue` statements
 * which may not be easy to abstract into a real function template.
 */
///@{

/** Macro for invoking a conventional separation routine.
 * In the body of Solver::cut_and_piv, this can be used to
 * invoke a straightforward separation routine which is just called and used
 * to update stats. See below for example usage.
 */
#define CUT_PIV_CALL(sep_ptr, find_fn, cuts_q, description)     \
try {                                                           \
    reset_separator(sep_ptr);                                   \
    if (call_separator([&sep_ptr, this]()                       \
                       { return sep_ptr->find_fn; },            \
                       sep_ptr->cuts_q(), piv, piv_stats)) {    \
        if (return_pivot(piv))                                  \
            return piv;                                         \
                                                                \
        if (restart_loop(piv, piv_stats.delta_ratio))           \
            continue;                                           \
    }                                                           \
} CMR_CATCH_PRINT_THROW(description, err);

/** Macro for invoking cut metamorphosis separation.
 * This can be used for a straightforward invocation of one of the
 * Sep::MetaCuts separation routines. See below for example usage.
 */
#define META_SEP_PIV_CALL(meta_type, description)                       \
try {                                                                   \
    reset_separator(meta_sep);                                          \
    if (call_separator([&meta_sep]()                                    \
                       { return meta_sep->find_cuts(meta_type); },      \
                       meta_sep->metacuts_q(), piv, piv_stats)) {       \
        if (return_pivot(piv))                                          \
            return piv;                                                 \
                                                                        \
        if (restart_loop(piv, piv_stats.delta_ratio))                   \
            continue;                                                   \
    }                                                                   \
}  CMR_CATCH_PRINT_THROW(description, err);

///@}

/**
 * @param do_price is edge pricing being performed for this instance. Used to
 * control Gomory cut separation.
 * @returns a LP::PivType corresponding to a non-degenerate pivot that was
 * obtainedfrom repeated rounds of cut generation and pivoting.
 */
PivType Solver::cut_and_piv(bool do_price)
{
    runtime_error err("Problem in Solver::cut_and_piv");

    PivType piv = PivType::Frac;
    int round = 0;
    bool &verbose = output_prefs.verbose;

    Data::SupportGroup &supp_data = core_lp.supp_data;

    unique_ptr<Sep::Separator> sep;
    unique_ptr<Sep::MetaCuts> meta_sep;

    if (output_prefs.prog_bar) {
        std::fill(p_bar.begin(), p_bar.end(), ' ');
        p_bar[0] = '[';
        p_bar[79] = ']';
    }

    while (true) {
        try { piv = core_lp.primal_pivot(); }
        CMR_CATCH_PRINT_THROW("initializing pivot and separator", err);

        if (return_pivot(piv))
            return piv;

        PivStats piv_stats(core_lp.get_objval());
        double &delta_ratio = piv_stats.delta_ratio;

        ++round;
        if (verbose) {
            printf("Round %d, first piv %.2f, ",
                   round, core_lp.get_objval());
            cout << piv << ", " << core_lp.num_rows() << " rows, "
                 << core_lp.num_cols() << " cols" << endl;
        }

        if (cut_sel.cutpool)
            CUT_PIV_CALL(sep, pool_sep(core_lp.ext_cuts), cutpool_q,
                         "doing pool sep");

        if (cut_sel.segment)
            CUT_PIV_CALL(sep, segment_sep(), segment_q, "doing segment sep");

        if (cut_sel.connect && !supp_data.connected) {
            if (verbose)
                cout << "\tAdding round of standard connect cuts...\n";
            int conrounds = 0;

            while (!supp_data.connected) {
                ++conrounds;
                try {
                    reset_separator(sep);
                    bool found_con =
                    call_separator([&sep]() { return sep->connect_sep(); },
                                   sep->connect_cuts_q(), piv, piv_stats);
                    if (!found_con)
                        throw logic_error("Disconnected w no connect cuts??");
                } CMR_CATCH_PRINT_THROW("doing connect cut loop", err);
            }

            if (return_pivot(piv)) {
                if (verbose)
                    cout << "\t...Tour pivot from "
                         << conrounds << " rounds of connect cuts" << endl;
                return piv;
            } else {
                if (verbose)
                    cout << "\t...Cut and pivoted for "
                         << conrounds << " rounds, solution now connected."
                         << endl;
                continue;
            }
        }

        if (cut_sel.fast2m)
            CUT_PIV_CALL(sep, fast2m_sep(), fastblossom_q, "doing fast2m sep");

        if (cut_sel.blkcomb)
            CUT_PIV_CALL(sep, blkcomb_sep(), blockcomb_q, "doing blkcomb sep");

        if (cut_sel.ex2m)
            CUT_PIV_CALL(sep, exact2m_sep(), exblossom_q, "doing exact 2m sep");

        if (cut_sel.simpleDP) {
            bool found_ex = false;
            int exrounds = 0;

            if (verbose)
                cout << "\tAdding round of standard exact SECs..." << endl;

            do {
                ++exrounds;
                try {
                    reset_separator(sep);
                    found_ex = call_separator([&sep]()
                                              { return sep->exsub_sep(); },
                                              sep->exact_sub_q(),
                                              piv, piv_stats);
                } CMR_CATCH_PRINT_THROW("doing exact subtour loop", err);
            } while (found_ex);

            if (return_pivot(piv)) {
                if (verbose)
                    cout << "\t...Tour pivot from " << exrounds
                         << " rounds of exact SEC generation" << endl;
                return piv;
            } else if (verbose)
                cout << "\t...Cut and pivoted for " << exrounds
                     << " rounds, CC says pivot now in subtour polytope..."
                     << endl;

            CUT_PIV_CALL(sep, simpleDP_sep(), simpleDP_q,
                         "doing simple DP sep");
        }

        using MetaType = Sep::MetaCuts::Type;

        if (cut_sel.tighten)
            CUT_PIV_CALL(meta_sep, tighten_cuts(), metacuts_q,
                         "tightening cuts");

        if (cut_sel.decker)
            META_SEP_PIV_CALL(MetaType::Decker, "doing Double Decker sep");

        if (cut_sel.handling)
            META_SEP_PIV_CALL(MetaType::Handling, "doing Handling sep");

        if (cut_sel.teething)
            META_SEP_PIV_CALL(MetaType::Teething, "doing Teething sep");

        if (cut_sel.tighten_pool)
            CUT_PIV_CALL(sep, tighten_pool(core_lp.ext_cuts), cutpool_q,
                         "tightening cut pool");

        if (cut_sel.consec1)
            CUT_PIV_CALL(sep, consec1_sep(core_lp.ext_cuts), consec1_q,
                         "doing consec1 comb sep");

        if (cut_sel.localcuts) {
            bool lc_restart = false;
            for (int chk = 8; chk <= Sep::LocalCuts::MaxChunkSize; ++chk) {
                reset_separator(sep);
                bool do_sphere = false;

                if (call_separator([&sep, chk, do_sphere]()
                                   { return sep->local_sep(chk, do_sphere); },
                                   sep->local_cuts_q(), piv, piv_stats)) {
                    if (return_pivot(piv))
                        return piv;

                    lc_restart = restart_loop(piv, delta_ratio);
                    if (lc_restart)
                        break;
                }
            }

            if (lc_restart)
                continue;

            for (int chk = 8; chk <= Sep::LocalCuts::MaxChunkSize; ++chk) {
                reset_separator(sep);
                bool do_sphere = true;

                if (call_separator([&sep, chk, do_sphere]()
                                   { return sep->local_sep(chk, do_sphere); },
                                   sep->local_cuts_q(), piv, piv_stats)) {
                    if (return_pivot(piv))
                        return piv;

                    lc_restart = restart_loop(piv, delta_ratio);
                    if (lc_restart)
                        break;
                }
            }

            if (lc_restart)
                continue;
        }

#if CMR_HAVE_SAFEGMI

        unique_ptr<Sep::SafeGomory> gmi_sep;

        if (cut_sel.safeGMI && !do_price)
            CUT_PIV_CALL(gmi_sep, find_cuts(), gomory_q, "doing safe GMI sep");

#endif
        double ph_init_prev = piv_stats.first_last_ratio;
        bool found_cuts = piv_stats.found_cuts;

        if (found_cuts && ph_init_prev >= 0.01)
            continue;

        if (output_prefs.prog_bar)
            cout << endl;

        if (verbose) {
            cout << "Tried all routines, found_cuts " << found_cuts
                 << ", returning " << piv << endl;
            piv_stats.report_extrema();
        }

        break;
    }

    return piv;
}

PivType Solver::abc_bcp(bool do_price)
{
    using BranchStat = ABC::BranchNode::Status;
    using ABC::BranchHistory;

    runtime_error err("Problem in Solver::abc_bcp");

    PivType piv = PivType::Frac;
    BranchHistory::iterator cur = branch_controller->next_prob();

    while (cur != branch_controller->get_history().end()) {
        cout << "\n";

        if (cur->stat == BranchStat::NeedsRecover) {
            cout << ABC::bnode_brief(*cur) << " needs feas recover"
                 << endl;
            if (do_price) {
                bool feasible = false;

                try { feasible = edge_pricer->feas_recover(); }
                CMR_CATCH_PRINT_THROW("doing feas_recover on node", err);

                if (feasible) {
                    cout << "Recovered infeasible BranchNode." << endl;
                    cur->stat = BranchStat::NeedsCut;
                } else {
                    cout << "Pricing pruned node by infeasibility." << endl;
                    cur->stat = BranchStat::Pruned;
                }
            } else {
                cout << "Infeasible sparse problem can be pruned" << endl;
                cur->stat = BranchStat::Pruned;
            }
        }

        if (cur->stat != BranchStat::Pruned)
            try {
                branch_controller->do_branch(*cur);
                if (core_lp.active_tourlen() < best_data.min_tour_value) {
                    cout << "Instated branch tour improves on best tour"
                         << endl;
                    core_lp.active_tour.best_update(best_data);
                    report_aug(Aug::Branch);
                    if (lb_fathom()) {
                        cout << "Branch tour is optimal by LB, "
                             << "terminating ABC search." << endl;
                        return PivType::FathomedTour;
                    }
                }
            } CMR_CATCH_PRINT_THROW("branching on current problem", err);

        if (cur->stat == BranchStat::NeedsPrice) {
            cout << "Price check: "
                 << ABC::bnode_brief(*cur) << " based on ";
            double opt_time = util::zeit();
            try {
                if (cur->price_basis) {
                    cout << "optimal estimate";
                    core_lp.copy_base(cur->price_basis->colstat,
                                      cur->price_basis->rowstat);
                } else {
                    cout << "high estimate";
                    core_lp.copy_start(core_lp.active_tour.edges());
                }
                cout << endl;

                core_lp.primal_opt();
            } CMR_CATCH_PRINT_THROW("optimizing for pricing", err);
            opt_time = util::zeit() - opt_time;

            double objval = core_lp.get_objval();

            if (output_prefs.verbose) {
                printf("\tPrimal opt objval %.2f in %.2fs", objval, opt_time);
                cout << endl;
            }

            cur->stat = BranchStat::NeedsCut;

            if (objval >= core_lp.global_ub() - 0.9) {
                if (!do_price) {
                    cout << "\tSparse problem prunable." << endl;
                    cur->stat = BranchStat::Pruned;
                } else {
                    Price::ScanStat pstat = Price::ScanStat::Full;
                    try {
                        pstat = edge_pricer->gen_edges(PivType::FathomedTour,
                                                       false);
                    } CMR_CATCH_PRINT_THROW("pricing edges", err);

                    if (pstat == Price::ScanStat::FullOpt) {
                        cur->stat = BranchStat::Pruned;
                        cout << "Problem pruned by LB, do not cut." << endl;
                    }
                }
            }
        }

        if (cur->stat == BranchStat::NeedsCut) {
            cout << "Cutting on " << ABC::bnode_brief(*cur) << endl;
            try {
                piv = cutting_loop(do_price, false, false);
            } CMR_CATCH_PRINT_THROW("cutting branch prob", err);

            if (piv == PivType::Frac)
                cur->stat = BranchStat::NeedsBranch;
            else if (piv == PivType::FathomedTour) {
                cur->stat = BranchStat::Pruned;
                cout << "\tPruned with opt objval "
                     << core_lp.get_objval() << endl;
                if (lb_fathom()) {
                    cout << "Terminating ABC search by lower bound." << endl;
                    return piv;
                }
            } else {
                cerr << "Pivot status " << piv << " in abc" << endl;
                throw err;
            }
        }

        if (cur->stat == BranchStat::NeedsBranch) {
            cout << "Branching on " << ABC::bnode_brief(*cur) << endl;
            try {
                branch_controller->split_prob(cur);
            } CMR_CATCH_PRINT_THROW("splitting branch problem", err);

            cur->stat = BranchStat::Done;
        }

        try { branch_controller->do_unbranch(*cur); }
        CMR_CATCH_PRINT_THROW("unbranching pruned problem", err);

        cur = branch_controller->next_prob();
    }

    return piv;
}

PivType Solver::frac_recover()
{
    runtime_error err("Problem in Solver::frac_recover");

    Data::SupportGroup &s_dat = core_lp.supp_data;
    int ncount = tsp_instance.node_count();
    vector<int> cyc;

    try { cyc.resize(ncount); } CMR_CATCH_PRINT_THROW("allocating cyc", err);

    double val = std::numeric_limits<double>::max();
    CCrandstate rstate;
    CCutil_sprand(tsp_instance.seed(), &rstate);

    if (CCtsp_x_greedy_tour_lk(tsp_instance.ptr(), ncount,
                               s_dat.support_ecap.size(),
                               &s_dat.support_elist[0],
                               &s_dat.support_ecap[0], &cyc[0], &val, true,
                               &rstate)) {
        cerr << "CCtsp_x_greedy_tour_lk failed.\n";
        throw err;
    }

    if (val >= core_lp.active_tourlen())
        return PivType::Frac;

    vector<Graph::Edge> new_edges;

    for (int i = 0; i < ncount; ++i) {
        EndPts e(cyc[i], cyc[(i + 1) % ncount]);
        if (core_graph.find_edge_ind(e.end[0], e.end[1]) == -1) {
            try {
                new_edges.emplace_back(e.end[0], e.end[1],
                                       tsp_instance.edgelen(e.end[0],
                                                            e.end[1]));
            } CMR_CATCH_PRINT_THROW("emplacing new edge", err);
        }
    }

    if (!new_edges.empty()) {
        try {
            core_lp.add_edges(new_edges, false);
        } CMR_CATCH_PRINT_THROW("adding edges not in tour", err);
    }

    if (output_prefs.verbose)
        cout << "Recovered to tour of length " << val << ", "
             << new_edges.size() << " new edges added" << endl;

    try { core_lp.set_active_tour(std::move(cyc)); }
    CMR_CATCH_PRINT_THROW("passing recover tour to core_lp", err);

    return LP::PivType::Tour;
}

void Solver::reset_separator(unique_ptr<Sep::Separator> &S)
{
    util::ptr_reset(S, core_graph.get_edges(), active_tour(),
                    core_lp.supp_data, karp_part,
                    tsp_instance.seed());
    S->filter_primal = !active_tour().tourless();
    S->verbose = output_prefs.verbose;
}

void Solver::reset_separator(std::unique_ptr<Sep::MetaCuts> &MS)
{
    util::ptr_reset(MS, core_lp.external_cuts(), graph_info().get_edges(),
                    active_tour(), core_lp.supp_data);
    MS->filter_primal = !active_tour().tourless();
    MS->verbose = output_prefs.verbose;
}

#if CMR_HAVE_SAFEGMI

void Solver::reset_separator(std::unique_ptr<Sep::SafeGomory> &GS)
{
    util::ptr_reset(GS, core_lp, active_tour().edges(),
                    core_lp.lp_edges);
    GS->filter_primal = !active_tour().tourless();
    GS->verbose = output_prefs.verbose;
}

#endif

}
