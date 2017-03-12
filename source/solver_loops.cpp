/**
 * @file
 * @brief Private Solver methods with recursive/loop calls.
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


/** Function template to pivot, find cuts, pivot back, add them, pivot again.
 * This function template is used to add one class of cuts at a time, attempt
 * to measure progress obtained, and return information on whether a further
 * separation routine should be called.
 * @tparam Qtype the queue representation of the cuts found by \p sepcall.
 * @param[in] sepcall a function which returns true if cuts are found, in which
 * case they are stored in \p sep_q.
 * @param[out] piv if cuts are found, this is the pivot type that occurred
 * after they were added.
 * @param[in,out] prev_val the value of the last pivot computed; will be set to
 * a new value if cuts are found.
 * @param[in,out] total_delta the sum of changes in objective values attained
 * by separation routines called thus far.
 * @param[out] delta_ratio if cuts are found, this is the ratio of the
 * difference between the new pivot value and \p prev_val divided by
 * \p tourlen.
 */
template<typename Qtype>
bool Solver::call_separator(const function<bool()> &sepcall, Qtype &sep_q,
                            PivType &piv, double &prev_val,
                            double &total_delta, double &delta_ratio,
                            double &lowest_piv)
{
    bool result = sepcall();
    if (result) {
        core_lp.pivot_back(true);
        core_lp.add_cuts(sep_q);
        piv = core_lp.primal_pivot();

        double new_val = core_lp.get_objval();
        double delta = std::abs(new_val - prev_val);
        double tourlen = core_lp.active_tourlen();

        total_delta += delta;
        if (new_val < lowest_piv)
            lowest_piv = new_val;

        double ph_delta = std::abs((new_val - prev_val)/(tourlen - prev_val));
        delta_ratio = ph_delta;

        if (output_prefs.verbose) {
            cout << "\t^^Cuts objval change " << prev_val << " -> "
                 << new_val << "\n";
            cout << "\tTotal delta " << total_delta << ", ratio "
                 << delta_ratio << ", pivot " << piv << endl;
        }
        prev_val = new_val;
    }

    return result;
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
    return piv == PivType::Tour || piv == PivType::FathomedTour;
}

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
    bool verbose = output_prefs.verbose;

    Data::SupportGroup &supp_data = core_lp.supp_data;

    unique_ptr<Sep::Separator> sep;
    unique_ptr<Sep::PoolCuts> pool_sep;


    while (true) {
        try {
            piv = core_lp.primal_pivot();
        } CMR_CATCH_PRINT_THROW("initializing pivot and separator", err);

        double total_delta = 0.0;
        double delta_ratio = 0.0;
        double prev_val = core_lp.get_objval();
        double lowest_piv = prev_val;

        if (return_pivot(piv))
            return piv;

        ++round;
        if (verbose)
            cout << "Cut and Piv round " << round << ", initial pivot value "
                 << prev_val << ", " << piv << endl;


        bool found_primal = false;

        if (cut_sel.cutpool && core_lp.external_cuts().pool_count() != 0)
            try {
                reset_separator(pool_sep);
                if (call_separator([&pool_sep]()
                                   { return pool_sep->find_cuts(); },
                                   pool_sep->pool_q(), piv,
                                   prev_val, total_delta, delta_ratio,
                                   lowest_piv)) {
                    found_primal = true;

                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("calling pool sep", err);

        bool found_seg = false;

        if (cut_sel.segment)
            try {
                reset_separator(sep);
                if (call_separator([&sep]() { return sep->segment_sep(); },
                                   sep->segment_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv)) {
                    found_primal = true;
                    found_seg = true;

                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("calling segment sep", err);

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
                                   sep->connect_cuts_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv);
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
            try {
                reset_separator(sep);
                if (call_separator([&sep]() { return sep->fast2m_sep(); },
                                   sep->fastblossom_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv)) {
                    found_primal = true;

                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("calling fast2m sep", err);

        if (cut_sel.blkcomb)
            try {
                reset_separator(sep);
                if (call_separator([&sep]() { return sep->blkcomb_sep(); },
                                   sep->blockcomb_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv)) {
                    found_primal = true;
                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("calling blkcomb sep", err);

        if (cut_sel.ex2m)
            try {
                reset_separator(sep);
                if (call_separator([&sep]() { return sep->exact2m_sep(); },
                                   sep->exblossom_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv)) {
                    found_primal = true;
                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("calling exact 2m sep", err);

        if (cut_sel.simpleDP && !found_seg && supp_data.connected)
            try {
                reset_separator(sep);
                if (call_separator([&sep]() { return sep->simpleDP_sep(); },
                                   sep->simpleDP_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv)) {
                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("calling simpleDP sep", err);

        if (verbose && cut_sel.safeGMI)
            cout << "\tBefore GMI, total_delta "
                 << total_delta << ", last ratio " << delta_ratio << endl;

        if (total_delta < 0.001 * (core_lp.active_tourlen() - lowest_piv)) {
            if (!found_primal && total_delta != 0.0)
                cout << "Total delta " << total_delta << ", active tourlen "
                     << core_lp.active_tourlen() << ", lowest piv "
                     << lowest_piv << ", setting found_primal false" << endl;
            found_primal = false;
        }

        if (found_primal)
            continue;

#if CMR_HAVE_SAFEGMI

        unique_ptr<Sep::SafeGomory> gmi_sep;

        if (cut_sel.safeGMI && !do_price)
            try {
                reset_separator(gmi_sep);
                if (call_separator([&gmi_sep]()
                                   { return gmi_sep->find_cuts(); },
                                   gmi_sep->gomory_q(), piv,
                                   prev_val, total_delta,
                                   delta_ratio, lowest_piv)) {
                    if (return_pivot(piv))
                        return piv;

                    if (restart_loop(piv, delta_ratio))
                        continue;
                }
            } CMR_CATCH_PRINT_THROW("doing safe GMI sep", err);

#endif

        if (verbose) {
            cout << "Tried all routines, returning " << piv
                 << " with total delta " << total_delta << endl;
        }

        break;
    }

    return piv;
}

PivType Solver::abc_dfs(int depth, bool do_price)
{
    using SplitIter = ABC::SplitIter;
    using NodeIter = ABC::BranchHistory::iterator;
    using BranchStat = ABC::BranchNode::Status;

    runtime_error err("Prolem in Solver::abc_dfs");

    PivType piv = PivType::Frac;

    try {
        piv = cutting_loop(do_price, false, false);
    } CMR_CATCH_PRINT_THROW("solving branch prob", err);

    if (piv != PivType::Frac) {
        if (piv == PivType::FathomedTour)
            return piv;
        else {
            cerr << "Pivot status " << piv << " in abc.\n";
            throw logic_error("Invalid pivtype in abc_dfs.");
        }
    }

    SplitIter branch_probs;

    try {
        branch_probs = dfs_brancher->next_level();
    } CMR_CATCH_PRINT_THROW("getting next level", err);

    for (NodeIter &P : branch_probs) {
        if (P == dfs_brancher->get_history().end())
            continue;

        try {
            dfs_brancher->do_branch(*P);
        } CMR_CATCH_PRINT_THROW("doing branch", err);

        bool call_again = true;

        if (P->stat == BranchStat::NeedsRecover) {
            cout << *P << " appears infeasible\n";
            cout << "\tActual verification should go here!!!" << endl;
            call_again = false;
            P->stat = BranchStat::Pruned;
        } else if (P->stat == BranchStat::NeedsPrice) {
            cout << *P << " looks optimal\n";
            if (!do_price) {
                cout << "Sparse edge set, problem can be pruned" << endl;
                call_again = false;
            } else {
                if (!P->price_basis)
                    throw runtime_error("No price basis for NeedsPrice node");
                try {
                    core_lp.copy_base(P->price_basis->colstat,
                                      P->price_basis->rowstat);
                    core_lp.primal_opt();
                    cout << "\tOpt objval : " << core_lp.get_objval() << endl;

                    Price::ScanStat price_stat =
                    edge_pricer->gen_edges(LP::PivType::FathomedTour, false);

                    if (price_stat == Price::ScanStat::FullOpt) {
                        cout << "Branch problem can be pruned." << endl;
                        P->stat = BranchStat::Pruned;
                        call_again = false;
                    } else {
                        core_lp.pivot_back(false);
                    }
                } CMR_CATCH_PRINT_THROW("processing NeedsPrice problem", err);
            }
        }

        if (call_again)
            piv = abc_dfs(++depth, do_price);

        try {
            dfs_brancher->do_unbranch(*P);
        } CMR_CATCH_PRINT_THROW("undoing branch", err);
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
                    core_lp.supp_data, karp_part);
    S->filter_primal = !branch_engaged;
    S->verbose = output_prefs.verbose;
}

void Solver::reset_separator(std::unique_ptr<Sep::PoolCuts> &PS)
{
    util::ptr_reset(PS, core_lp.ext_cuts,
                    core_graph.get_edges(), active_tour().edges(),
                    core_lp.supp_data);
    PS->filter_primal = !branch_engaged;
    PS->verbose = output_prefs.verbose;
}

#if CMR_HAVE_SAFEGMI

void Solver::reset_separator(std::unique_ptr<Sep::SafeGomory> &GS)
{
    util::ptr_reset(GS, core_lp, active_tour().edges(),
                    core_lp.lp_edges);
    GS->filter_primal = !branch_engaged;
    GS->verbose = output_prefs.verbose;
}

#endif

}
