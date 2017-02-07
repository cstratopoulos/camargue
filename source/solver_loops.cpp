/**
 * @file
 * @brief Implementation of private Solver methods with recursive calls.
 */

#include "config.hpp"
#include "solver.hpp"
#include "separator.hpp"
#include "brancher.hpp"

#if CMR_HAVE_SAFEGMI
#include "safeGMI.hpp"
#endif

#include "err_util.hpp"

#include <array>
#include <iostream>
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
 * @param[in] core_lp the LP relaxation for performing pivots and adding cuts.
 * @param[in] tourlen the length of the best known tour.
 * @param[in/out] prev_val the value of the last pivot computed; will be set to
 * a new value if cuts are found.
 * @param[in/out] total_delta the sum of changes in objective values attained
 * by separation routines called thus far.
 * @param[out] delta_ratio if cuts are found, this is the ratio of the 
 * difference between the new pivot value and \p prev_val divided by 
 * \p tourlen.
 * @param[out] num_pruned number of cuts pruned from the LP after a pivot call.
 * @returns the value of \p sepcall, i.e., true iff cuts are found.
 */
template<class Qtype>
bool call_separator(const function<bool()> &sepcall, const Qtype &sep_q,
                    PivType &piv, LP::CoreLP &core_lp,
                    const double tourlen, double &prev_val,
                    double &total_delta, double &delta_ratio, int &num_pruned)
{
    bool result = sepcall();
    int num_rows = core_lp.num_rows();
    if (result) {
        core_lp.pivot_back();
        core_lp.add_cuts(sep_q);
        piv = core_lp.primal_pivot();

        if (piv == PivType::Tour) {
            num_pruned = num_rows - core_lp.num_rows();
        }

        double new_val = core_lp.get_objval();
        double delta = std::abs(new_val - prev_val);
        prev_val = new_val;
        
        total_delta += delta;
        delta_ratio = (delta / tourlen);
    }

    return result;
}

PivType Solver::cut_and_piv(int &round, int &num_pruned, bool do_price)
{
    runtime_error err("Problem in Solver::cut_and_piv");
    bool silent = true;
    
    const double tourlen = core_lp.get_objval();
    double prev_val = tourlen;
    double total_delta = 0.0;
    double delta_ratio = 0.0;
    
    PivType piv;

    Data::SupportGroup &supp_data = core_lp.supp_data;
    unique_ptr<Sep::Separator> sep;

    ++round;
    if (!silent)
        cout << "\nRound " << round << "\n";

    int num_rows = core_lp.num_rows();

    try {
        piv = core_lp.primal_pivot();
        sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                supp_data, karp_part, TG);
    } CMR_CATCH_PRINT_THROW("initializing pivot and separator", err);

    if (piv == PivType::Tour || piv == PivType::FathomedTour) {
        num_pruned = num_rows - core_lp.num_rows();
        return piv;
    }

    bool found_seg = false;
    bool found_primal = false;

    if (cut_sel.segment)
        try {
            if (call_separator([&sep]() { return sep->segment_sep(); },
                               sep->segment_q(), piv, core_lp,
                               tourlen, prev_val, total_delta, delta_ratio,
                               num_pruned)) {
                found_primal = true;
                found_seg = true;
                if (piv == PivType::Tour || piv == PivType::FathomedTour)
                    return piv;

                if (piv == PivType::Subtour || delta_ratio > Eps::SepRound)
                    return cut_and_piv(round, num_pruned, do_price);

                sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                        supp_data, karp_part,
                                                        TG);
            }
        } CMR_CATCH_PRINT_THROW("calling segment sep", err);

    if (cut_sel.fast2m)
        try {
            if (call_separator([&sep]() { return sep->fast2m_sep(); },
                               sep->fastblossom_q(), piv, core_lp,
                               tourlen, prev_val, total_delta, delta_ratio,
                               num_pruned)) {
                found_primal = true;
                if (piv == PivType::Tour || piv == PivType::FathomedTour)
                    return piv;

                if (piv == PivType::Subtour || delta_ratio > Eps::SepRound)
                    return cut_and_piv(round, num_pruned, do_price);

                sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                        supp_data, karp_part,
                                                        TG);
            }
        } CMR_CATCH_PRINT_THROW("calling fast2m sep", err);

    if (cut_sel.blkcomb)
        try {
            if (call_separator([&sep]() { return sep->blkcomb_sep(); },
                               sep->blockcomb_q(), piv, core_lp,
                               tourlen, prev_val, total_delta, delta_ratio,
                               num_pruned)) {
                found_primal = true;
                if (piv == PivType::Tour || piv == PivType::FathomedTour)
                    return piv;

                if (piv == PivType::Subtour || found_seg ||
                    !supp_data.connected)
                    return cut_and_piv(round, num_pruned,  do_price);

                sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                        supp_data, karp_part,
                                                        TG);
            }
        } CMR_CATCH_PRINT_THROW("calling blkcomb sep", err);

    if (cut_sel.simpleDP)
        try {
            if (!found_seg && supp_data.connected &&
                call_separator([&sep]() { return sep->simpleDP_sep(); },
                               sep->simpleDP_q(), piv, core_lp,
                               tourlen, prev_val, total_delta, delta_ratio,
                               num_pruned)) {
                if (piv == PivType::Tour || piv == PivType::FathomedTour)
                    return piv;

                if (total_delta >= Eps::Zero)
                    return cut_and_piv(round, num_pruned,  do_price);
            }
        } CMR_CATCH_PRINT_THROW("calling simpleDP sep", err);

    if (cut_sel.connect) {
        if (!found_primal && !supp_data.connected) {
            int num_add = 0;
            while (!supp_data.connected) {
                try {
                    if (call_separator([&sep]() { return sep->connect_sep(); },
                                       sep->connect_cuts_q(), piv, core_lp,
                                       tourlen, prev_val, total_delta,
                                       delta_ratio, num_pruned)) {
                        num_add += sep->connect_cuts_q().size();
                        sep = util::make_unique<Sep::Separator>(graph_data,
                                                                best_data,
                                                                supp_data,
                                                                karp_part, TG);

                    } else {
                        throw logic_error("Disconnected w no connect cuts??");
                    }                
                } CMR_CATCH_PRINT_THROW("doing connect cut loop", err);
            }
        
            if (piv == PivType::Tour || piv == PivType::FathomedTour) {
                return piv;
            } else {
                return cut_and_piv(round, num_pruned,  do_price);
            }
        } else if (!silent) {
            cout << "\tcuts: " << found_primal << ",connected "
                 << supp_data.connected << "\n";
        }
    }

    if (total_delta < Eps::Zero)
        found_primal = false;

    if (found_primal) {
        return cut_and_piv(round, num_pruned, do_price);
    }

#if CMR_HAVE_SAFEGMI
    
    if (cut_sel.safeGMI)
        if (!do_price) {
            try {
                vector<double> lp_x = core_lp.lp_vec();
                Sep::SafeGomory gmi_sep(core_lp,
                                        core_lp.tour_base.best_tour_edges,
                                        lp_x);
                
                if (call_separator([&gmi_sep]() { return gmi_sep.find_cuts(); },
                                   gmi_sep.gomory_q(), piv, core_lp,
                                   tourlen, prev_val, total_delta,
                                   delta_ratio, num_pruned)) {
                    if (piv == PivType::Tour || piv == PivType::FathomedTour)
                        return piv;

                    if (total_delta > Eps::Zero)
                        return cut_and_piv(round, num_pruned, do_price);
                }
            } CMR_CATCH_PRINT_THROW("doing safe GMI sep", err);
        }
    
#endif


    if (!silent)
        cout << "\tTried all routines, returning " << piv << "\n";
    return piv;    
}

PivType Solver::abc_dfs(int depth, bool do_price)
{
    using Problem = ABC::Problem;
    using Ptype = Problem::Type;
    using ProbArray = std::array<Problem, 2>;
    
    runtime_error err("Prolem in Solver::abc_dfs");

    PivType piv = PivType::Frac;

    if (depth > 0)
        try {
            piv = cutting_loop(do_price, false);
        } CMR_CATCH_PRINT_THROW("solving branch prob", err);

    if (piv != PivType::Frac) {
        if (piv == PivType::FathomedTour)
            return piv;
        else {
            cerr << "Pivot status " << piv << " in abc.\n";
            throw logic_error("Invalid pivtype in abc_dfs.");
        }
    }

    ProbArray branch_probs = brancher->next_level();
    const double &tourlen = best_data.min_tour_value;
    const vector<double> &tour_vec = tour_basis().best_tour_edges;


    for (Problem &P : branch_probs) {
        try {
            cout << "\n\tSearch depth " << depth << "\n";
            brancher->do_branch(P);
            if (P.type == Ptype::Affirm)
                core_lp.copy_start(tour_vec);
            else
                core_lp.copy_base(*P.contra_base);

            core_lp.factor_basis();
        }
        CMR_CATCH_PRINT_THROW("doing branch", err);

        double score = P.rank.first;
        double estimate = P.rank.second;
        bool call_again = true;

        if (score == 2) {
            cout << "\tProblem appears infeasible.\n";
            cout << "\tVerification should go here!!!\n";
            call_again = false;
        } else if (score == 1) {
            if (estimate >= tourlen - 0.9) {
                cout << "\tProblem appears prunable.\n";
                if (!do_price) {
                    cout << "\tSparse problem/no price, prune search.\n";
                    call_again = false;
                }
            }
        }

        if (call_again)
            piv = abc_dfs(++depth, do_price);        

        try {
            cout << "\n\tSearch depth " << depth << "\n";
            brancher->undo_branch(P);
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

    double val = DoubleMax;
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

    if (val >= best_data.min_tour_value) {
        //cout << "\tUnable to augment from frac vector.\n";
        return PivType::Frac;
    } else
        cout << "\tFound improved tour from frac vector:\n"
             << "\t\t" << best_data.min_tour_value << " --> " << val << "\n";
        
    vector<Graph::Edge> new_edges;
    Graph::CoreGraph &graph = graph_data.core_graph;

    for (int i = 0; i < ncount; ++i) {
        EndPts e(cyc[i], cyc[(i + 1) % ncount]);
        if (graph.find_edge_ind(e.end[0], e.end[1]) == -1) {
            try {
                new_edges.emplace_back(e.end[0], e.end[1],
                                       tsp_instance.edgelen(e.end[0],
                                                            e.end[1]));
            } CMR_CATCH_PRINT_THROW("emplacing new edge", err);
        }
    }

    if (!new_edges.empty()) {
        int orig_rowcount = core_lp.num_rows();
        try {
            if (cut_sel.safeGMI)
                core_lp.purge_gmi();
            core_lp.add_edges(new_edges);
        } CMR_CATCH_PRINT_THROW("adding edges not in tour", err);
        int new_rowcount = core_lp.num_rows();
        cout << "\tRecover tour contains " << new_edges.size() << " new edges, "
             << (orig_rowcount - new_rowcount) << " gmi cuts purged.\n";
    }

    vector<double> &lp_edges = core_lp.lp_edges;

    for (int i = 0; i < lp_edges.size(); ++i)
        lp_edges[i] = 0.0;

    //prepping for basis rebuild/handle aug
    for (int i = 0; i < ncount; ++i) {
        EndPts e(cyc[i], cyc[(i + 1) % ncount]);
        int ind = graph.find_edge_ind(e.end[0], e.end[1]);
        if (ind == -1) {
            cerr << "Tour edge " << e.end[0] << ", " << e.end[1]
                 << " still not in graph\n";
            throw err;
        }
        lp_edges[ind] = 1.0;
    }

    //used as the tour nodes vector in handle_aug
    graph_data.island = std::move(cyc);

    try {
        core_lp.copy_start(lp_edges);
        core_lp.factor_basis();
        core_lp.handle_aug(); //instates the tour stored in lp_edges
    } CMR_CATCH_PRINT_THROW("rebuilding/augmenting for x-tour", err);
    
    return LP::PivType::Tour;
}

}
