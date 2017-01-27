#include "solver.hpp"
#include "separator.hpp"
#include "timer.hpp"
#include "err_util.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>

using std::abs;

using std::cout;
using std::cerr;
using std::endl;

using std::string;

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::unique_ptr;

using std::vector;


namespace CMR {

using CutType = Sep::HyperGraph::Type;
using PivType = LP::PivType;
namespace Eps = Epsilon;

inline static int make_seed(const int seed)
{
    return (seed > 0) ? seed : (int) util::real_zeit();
}

Solver::Solver(const string &fname, const int seed, const OutPrefs outprefs)
try : tsp_instance(fname, make_seed(seed)),
      karp_part(tsp_instance),
      graph_data(tsp_instance), best_data(tsp_instance, graph_data),
      core_lp(graph_data, best_data),
      output_prefs(outprefs)
    
{
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver TSPLIB constructor failed.");
}

Solver::Solver(const string &fname, const string &tour_fname,
               const int seed, const OutPrefs outprefs)
try : tsp_instance(fname, make_seed(seed)),
      karp_part(tsp_instance),
      graph_data(tsp_instance),
      best_data(tsp_instance, graph_data, tour_fname),
      core_lp(graph_data, best_data),
      output_prefs(outprefs)
{
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver TSPLIB/tour constructor failed.");
}

Solver::Solver(const int seed, const int node_count,
               const int gridsize, const OutPrefs outprefs)
try : tsp_instance(make_seed(seed), node_count, gridsize),
      karp_part(tsp_instance),
      graph_data(tsp_instance), best_data(tsp_instance, graph_data),
      core_lp(graph_data, best_data),
      output_prefs(outprefs)
{
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver random constructor failed.");
}

void Solver::report_piv(PivType piv, int round, bool full_opt)
{
    cout << "\n\t\tRound " << round << "\n"
         << "\tLP objective value: " << core_lp.get_objval()
         << ", dual feasible: " << core_lp.dual_feas() << "\n";
    cout << "\t" << core_lp.num_rows() << " rows, "
         << core_lp.num_cols() << " cols in LP.\n";
    
    switch (piv) {
    case PivType::FathomedTour:
        if (full_opt) {
            cout << "\tCurrent tour is optimal\n"
                 << "\t*************************\n";
        } else {
            cout << "\tTour optimal for edge set, pricing edges...\n";
        }
        break;
    case PivType::Tour:
        cout << "\tAugmented to new tour.\n";
        break;
    default:
        cout << "\t Pivot status: \t" << piv << "\n";
    }
    
    cout << endl;
}

PivType Solver::cutting_loop(bool do_price)
{
    runtime_error err("Problem in Solver::cutting_loop.");

    if (do_price)
        try {
            edge_pricer = util::make_unique<Price::Pricer>(core_lp,
                                                           tsp_instance,
                                                           graph_data);
        } CMR_CATCH_PRINT_THROW("instantiating/allocating Pricer", err);
    
    PivType piv = PivType::Frac;
    int round = 0;
    int auground = 0;
    int stag_round = 0;

    vector<int> &tour_edges = best_data.best_tour_edges;
    const vector<Graph::Edge> &edges = graph_data.core_graph.get_edges();
    vector<int> &perm = best_data.perm;

    try {
        TG = Graph::TourGraph(tour_edges, edges, perm);        
    } CMR_CATCH_PRINT_THROW("allocating tour graph", err);


    CMR::Timer timer(tsp_instance.problem_name());
    
    timer.start();

    

    while (true) {
        ++auground;

        try {
            piv = cut_and_piv(round, stag_round);
        } CMR_CATCH_PRINT_THROW("invoking cut and piv", err);

        if (piv == PivType::FathomedTour) {
            report_piv(piv, round, !do_price);
            
            if (do_price) {
                try {
                    if (edge_pricer->gen_edges(piv) == Price::ScanStat::Full) {
                        core_lp.rebuild_basis();
                        core_lp.pivot_back();
                        continue;
                    } else
                        break;
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }            
            break;
        }

        if (piv == PivType::Tour) {
            report_piv(piv, round, false);

            if (do_price) {
                try {
                    edge_pricer->gen_edges(piv);
                    core_lp.rebuild_basis();
                    core_lp.pivot_back();
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }
            
            try {
                TG = Graph::TourGraph(tour_edges, edges, perm);
            } CMR_CATCH_PRINT_THROW("updating tour graph", err);
            
            continue;
        }

        report_piv(piv, round, false);
        cout << "\tNo cuts found.\n";
        break;
    }

    timer.stop();
    timer.report(true);

    if (piv == PivType::FathomedTour && do_price)
        report_piv(piv, round, true);

    int subcount = 0;
    int combcount = 0;
    int dpcount = 0;
    
    for (const Sep::HyperGraph &H : core_lp.ext_cuts.get_cuts()) {
        CutType t = H.cut_type();
        if (t == CutType::Subtour)
            ++subcount;
        else if (t == CutType::Comb)
            ++combcount;
        else if (t == CutType::Domino)
            ++dpcount;
    }

    cout << "\t" << subcount << " SECs, " << combcount << " combs/blossoms, "
         << dpcount << " dp cuts.\n";
    cout << "\n";
    return piv;
}

PivType Solver::abc(bool do_price)
{
    runtime_error err("Problem in Solver::abc.");

    PivType piv = cutting_loop(do_price);
    
    if (piv != PivType::Frac) {
        if (piv == PivType::FathomedTour)
            return piv;
        else {
            cerr << "Pivot status " << piv << " in abc.\n";
            throw logic_error("Invalid pivot type for running Solver::abc.");
        }            
    }

    try {
        brancher = util::make_unique<ABC::Brancher>(core_lp,
                                                    graph_data.core_graph
                                                    .get_edges(), tour_basis(),
                                                    best_data.min_tour_value,
                                                    ABC::ContraStrat::Fix);
    } CMR_CATCH_PRINT_THROW("allocating/instantiating Brancher", err);

    using ProbStat = ABC::Problem::Status;
    ABC::Problem prob = brancher->next_prob();

    cout << "\n\n\t\t///Beginning ABC search\n\n";

    while (!brancher->solved(prob)) {
        cout << "\tBRANCH TASK: " << prob << "\n";
        piv = cutting_loop(do_price);

        cout << "\tTASK STATUS: ";
        if (piv == PivType::FathomedTour) {
            brancher->pop_problem(ProbStat::Pruned);
        } else if (piv == PivType::Frac) {
            brancher->pop_problem(ProbStat::Seen);
        } else {
            cerr << "Pivot status " << piv << " in abc loop.\n";
            throw runtime_error("Invalid piv stat for further branching.");
        }
        
        cout << "Calling next prob...." << endl;
        prob = brancher->next_prob();
    }

    cout << "\tABC search complete.\n";
    
    return piv;    
}

PivType Solver::cut_and_piv(int &round, int &stag_rounds)
{
    runtime_error err("Problem in Solver::cut_and_piv");
    bool silent = true;
    
    const double tourlen = core_lp.get_objval();
    double prev_val = tourlen;
    double total_delta = 0.0;
    
    PivType piv;

    Data::SupportGroup &supp_data = core_lp.supp_data;
    unique_ptr<Sep::Separator> sep;

    ++round;
    if (!silent)
        cout << "\nRound " << round << "\n";
    if (stag_rounds > 10) {
        if (!silent)
            cout << "\tTerminating due to stagnation.\n";
        return PivType::Frac;
    }

    try {
        piv = core_lp.primal_pivot();
        sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                supp_data, karp_part, TG);
    } CMR_CATCH_PRINT_THROW("initializing pivot and separator", err);

    if (piv == PivType::Tour || piv == PivType::FathomedTour)
        return piv;

    bool found_seg = false;
    bool found_primal = false;

    try {
        if (sep->segment_sep()) {
            found_seg = true;
            found_primal = true;
        
            core_lp.pivot_back();
            core_lp.add_cuts(sep->segment_q());
            piv = core_lp.primal_pivot();

            double new_val = core_lp.get_objval();
            double delta = abs(new_val - prev_val);
            total_delta += delta;

            if (!silent)
                cout << "\tAdded " << sep->segment_q().size() << " segments, "
                     << piv << ", delta ratio " << (delta / tourlen) << "\n";

            if (piv == PivType::Tour || piv == PivType::FathomedTour)
                return piv;

            if (piv == PivType::Subtour || (delta / tourlen > Eps::SepRound))
                return cut_and_piv(round, stag_rounds);

            prev_val = new_val;
            sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                    supp_data, karp_part, TG);
        } else if (!silent) {
            cout << "\tNo segments.\n";
        }
    } CMR_CATCH_PRINT_THROW("doing segment sep", err);
    
    try {
        if (sep->fast2m_sep()) {
            found_primal = true;
            core_lp.pivot_back();
            core_lp.add_cuts(sep->fastblossom_q());
            piv = core_lp.primal_pivot();

            double new_val = core_lp.get_objval();
            double delta = abs(new_val - prev_val);
            total_delta += delta;

            if (!silent)
                cout << "\tAdded " << sep->fastblossom_q().size()
                     << " blossoms, " << piv << ", delta ratio "
                     << (delta / tourlen) << "\n";

            if (piv == PivType::Tour || piv == PivType::FathomedTour)
                return piv;

            if (piv == PivType::Subtour || (delta / tourlen > Eps::SepRound))
                return cut_and_piv(round, stag_rounds);

            prev_val = new_val;
            sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                    supp_data, karp_part, TG);
        } else if (!silent) {
            cout << "\tNo fastblossoms\n";
        }
    } CMR_CATCH_PRINT_THROW("doing fastblossom sep", err);

    try {
        if (sep->blkcomb_sep()) {
            found_primal = true;
            core_lp.pivot_back();
            core_lp.add_cuts(sep->blockcomb_q());
            piv = core_lp.primal_pivot();

            double new_val = core_lp.get_objval();
            double delta = abs(new_val - prev_val);
            total_delta += delta;

            if (!silent)
                cout << "\tAdded " << sep->blockcomb_q().size()
                     << " blk combs, " << piv << ", delta ratio "
                     << (delta / tourlen) << "\n";

            if (piv == PivType::Tour || piv == PivType::FathomedTour)
                return piv;

            if (piv == PivType::Subtour || found_seg || !supp_data.connected)
                return cut_and_piv(round, stag_rounds);

            prev_val = new_val;
            sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                    supp_data, karp_part, TG);
        } else if (!silent) {
            cout << "\tno blkcomb\n";
        }
    } CMR_CATCH_PRINT_THROW("doing block comb sep", err);

    try {
        if (!found_seg && supp_data.connected && sep->simpleDP_sep()) {
            found_primal = true;
            core_lp.pivot_back();
            core_lp.add_cuts(sep->simpleDP_q());
            piv = core_lp.primal_pivot();

            double new_val = core_lp.get_objval();
            double delta = abs(new_val - prev_val);
            total_delta += delta;

            if (!silent)
                cout << "\tAdded " << sep->simpleDP_q().size() << " dominos, "
                     << piv << ", delta ratio " << (delta / tourlen) << "\n";

            if (piv == PivType::Tour || piv == PivType::FathomedTour)
                return piv;
            else
                return cut_and_piv(round, stag_rounds);

            prev_val = new_val;
            sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                    supp_data, karp_part,
                                                    TG);
        } else if (!silent) {
            cout << "\tNo dominos.\n";
        }
    } CMR_CATCH_PRINT_THROW("doing simple DP sep", err);

    if (!found_primal && !supp_data.connected) {
        int num_add = 0;
        while (!supp_data.connected) {
            try {
                core_lp.pivot_back();
            
                if (!sep->connect_sep())
                    throw logic_error("No connect cuts in subtour solution.");
                total_delta += 1.0;

                core_lp.add_cuts(sep->connect_cuts_q());
                
                num_add += sep->connect_cuts_q().size();
                
                piv = core_lp.primal_pivot();
                prev_val = core_lp.get_objval();
                sep = util::make_unique<Sep::Separator>(graph_data, best_data,
                                                        supp_data, karp_part,
                                                        TG);
            } CMR_CATCH_PRINT_THROW("doing connect cut loop", err);
        }
        if (!silent)
            cout << "\tAdded " << num_add
                 << " connect cuts, piv " << piv << "\n";
        if (piv == PivType::Tour || piv == PivType::FathomedTour)
            return piv;
        else
            return cut_and_piv(round, stag_rounds);
    } else if (!silent) {
        cout << "\tcuts: " << found_primal << ",connected "
             << supp_data.connected << "\n";
    }

    if (found_primal) {
        if (total_delta < Eps::Zero)
            ++stag_rounds;
        else
            stag_rounds = 0;
        return cut_and_piv(round, stag_rounds);
    }

    if (!silent)
        cout << "\tTried all routines, returning " << piv << "\n";
    return piv;    
}

}
