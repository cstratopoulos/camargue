/**
 * @file
 * @brief Implementation of higher level or public Solver methods.
 */

#include "config.hpp"
#include "solver.hpp"

#include "timer.hpp"
#include "err_util.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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
{} catch (const exception &e) {
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
{} catch (const exception &e) {
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
{} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver random constructor failed.");
}

void Solver::report_lp(PivType piv)
{
    int rowcount = core_lp.num_rows();
    cout << "\tPivot status:\t" << piv << "\n"
         << "\tLP objective value: " << core_lp.get_objval()
         << ", dual feasible: " << core_lp.dual_feas() << "\n";
    cout << "\t" << rowcount << " rows, "
         << core_lp.num_cols() << " cols in LP.\n\n";
}

/**
 * Prints the progress in tour improvement.
 * @param piv_aug should be true if the tour was obtained from a non-degenerate
 * primal pivot, false if it was obtained from frac_recover.
 */
void Solver::report_aug(bool piv_aug)
{
    cout << "\tTour " << ++num_augs << ": "
         << core_lp.get_objval() << ", augmented from "
         << (piv_aug ? "primal pivot" : "x-heuristic") << "\n";
}

void Solver::report_cuts()
{
    int subcount = 0;
    int combcount = 0;
    int dpcount = 0;
    int gmicount = 0;
    
    for (const Sep::HyperGraph &H : core_lp.ext_cuts.get_cuts()) {
        CutType t = H.cut_type();
        if (t == CutType::Subtour)
            ++subcount;
        else if (t == CutType::Comb)
            ++combcount;
        else if (t == CutType::Domino)
            ++dpcount;
        else if (t == CutType::Non)
            ++gmicount;
    }

    cout << "\t" << subcount << " SECs, " << combcount
         << " combs/blossoms, " << dpcount << " dp cuts, "
         << gmicount << " GMI cuts. \n\t("
         << (core_lp.num_rows() - tsp_instance.node_count())
         << " cuts total, " << core_lp.num_cols() << " cols).\n";
    cout << "\n";
}


PivType Solver::cutting_loop(bool do_price, bool try_recover, bool pure_cut)
{
    runtime_error err("Problem in Solver::cutting_loop");

    if (do_price)
        try {
            edge_pricer = util::make_unique<Price::Pricer>(core_lp,
                                                           tsp_instance,
                                                           graph_data);
        } CMR_CATCH_PRINT_THROW("instantiating/allocating Pricer", err);
    
    PivType piv = PivType::Frac;
    int round = 0;
    int auground = 0;

    vector<int> &tour_edges = best_data.best_tour_edges;
    const vector<Graph::Edge> &edges = graph_data.core_graph.get_edges();
    vector<int> &perm = best_data.perm;

    try {
        TG = Graph::TourGraph(tour_edges, edges, perm);        
    } CMR_CATCH_PRINT_THROW("allocating tour graph", err);


    CMR::Timer timer(tsp_instance.problem_name() + " pure cut");
    
    timer.start();    

    while (true) {
        ++auground;
        int num_pruned = 0;

        try {
            piv = cut_and_piv(round, num_pruned, do_price);
        } CMR_CATCH_PRINT_THROW("invoking cut and piv", err);

        if (piv == PivType::FathomedTour) {
            if (do_price) {
                cout << "\tTour optimal for edge set...";
                try {
                    if (edge_pricer->gen_edges(piv) == Price::ScanStat::Full) {
                        core_lp.factor_basis();
                        //core_lp.rebuild_basis();
                        continue;
                    } else
                        break;
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }            
            break;
        }

        if (piv == PivType::Tour) {
            report_aug(true);

            if (do_price) {
                try {
                    edge_pricer->gen_edges(piv);
                    core_lp.factor_basis();
                    //core_lp.rebuild_basis();
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }
            
            try {
                TG = Graph::TourGraph(tour_edges, edges, perm);
            } CMR_CATCH_PRINT_THROW("updating tour graph", err);
            
            continue;
        }

        if (try_recover) {
            try {
                int prev_rows = core_lp.num_rows();
                piv = frac_recover();
                num_pruned = prev_rows - core_lp.num_rows();
                if (piv == PivType::Tour) {
                    report_aug(false);
                    TG = Graph::TourGraph(tour_edges, edges, perm);

                    continue;
                }
            } CMR_CATCH_PRINT_THROW("trying to recover from frac tour", err);
        }

        if (pure_cut) {
            cout << "\tNo cuts found.\n";
        }
        break;
    }

    timer.stop();

    if (pure_cut) {
        if (piv == PivType::FathomedTour)
            report_lp(piv);
        report_cuts();

        timer.report(true);        
    }
    return piv;
}

PivType Solver::abc(bool do_price)
{
    runtime_error err("Problem in Solver::abc");

    PivType piv = PivType::Frac;

    try { piv = cutting_loop(do_price, true, true); }
    CMR_CATCH_PRINT_THROW("running cutting_loop", err);
    
    if (piv != PivType::Frac) {
        if (piv == PivType::FathomedTour) {            
            return piv;
        }
        else {
            cerr << "Pivot status " << piv << " in abc.\n";
            throw logic_error("Invalid pivot type for running Solver::abc.");
        }            
    }

    cout << "\tCommencing ABC search....\n";
    Timer abct(tsp_instance.problem_name() + " ABC search");
    abct.start();

    try {
        brancher = util::make_unique<ABC::Brancher>(core_lp,
                                                    graph_data.core_graph
                                                    .get_edges(), tour_basis(),
                                                    best_data.min_tour_value,
                                                    ABC::ContraStrat::Fix);
    } CMR_CATCH_PRINT_THROW("allocating/instantiating Brancher", err);

    try { piv = abc_dfs(0, do_price); }
    CMR_CATCH_PRINT_THROW("running abc_dfs", err);

    abct.stop();

    cout << "\n\tABC search completed, optimal tour has length "
         << best_data.min_tour_value << "\n";
    report_cuts();

    abct.report(true);
    
    
    return piv;    
}

}
