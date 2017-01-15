#include "solver.hpp"
#include "separator.hpp"
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

using std::vector;


namespace CMR {

using CutType = Sep::HyperGraph::Type;

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

void Solver::report_piv(LP::PivType piv, int round)
{
    cout << "\n\t\tRound " << round << endl;
    
    switch (piv) {
    case LP::PivType::FathomedTour:
        cout << "\tTour optimal for edge set\n"
             << "\t*************************\n"
             << "\tLP optimal obj val: " << core_lp.get_objval()
             << ", dual feas: " << core_lp.dual_feas() << "\n";
        break;
    case LP::PivType::Tour:
        cout << "\tAugmented to tour of length "
             << core_lp.get_objval() << "\n";
        break;
    default:
        cout << "\t Pivot status: \t" << LP::piv_string(piv) << "\n"
             << "\t obj val: \t"
             << core_lp.get_objval() << ", dual feas: "
             << core_lp.dual_feas() << "\n";
    }

    cout << "\t" << core_lp.num_rows() << " rows, "
         << core_lp.num_cols() << " cols in LP.\n";
    
    cout << endl;
}

LP::PivType Solver::cutting_loop(bool do_price)
{
    runtime_error err("Problem in Solver::cutting_loop.");

    if (do_price)
        try {
            edge_pricer = util::make_unique<Price::Pricer>(core_lp,
                                                           tsp_instance,
                                                           graph_data);
        } CMR_CATCH_PRINT_THROW("Couldn't instantiate Pricer", err);
    
    LP::PivType piv = LP::PivType::Frac;
    int round = 0;
    int auground = 0;

    vector<int> &tour_edges = best_data.best_tour_edges;
    const vector<Graph::Edge> &edges = graph_data.core_graph.get_edges();
    vector<int> &perm = best_data.perm;

    Graph::TourGraph TG;

    try {
        TG = Graph::TourGraph(tour_edges, edges, perm);        
    } CMR_CATCH_PRINT_THROW("allocating tour graph", err);


    CMR::Timer timer(tsp_instance.problem_name());
    
    timer.start();

    while (true) {
        ++round;
        ++auground;
        int rowcount = core_lp.num_rows();

        try {
            piv = core_lp.primal_pivot();
        } CMR_CATCH_PRINT_THROW("performing primal pivot", err);
        
        if (piv == LP::PivType::FathomedTour) {
            report_piv(piv, round);
            
            if (do_price) {
                try {
                    if (edge_pricer->gen_edges(piv) == Price::ScanStat::Full){
                        continue;
                    } else
                        break;
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }
            
            break;
        }

        if (piv == LP::PivType::Tour) {
            report_piv(piv, round);
            cout << "\tPruned " << (rowcount - core_lp.num_rows())
                 << " rows from the LP." << endl;

            if (do_price)
                try {
                    edge_pricer->gen_edges(piv);
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            
            piv = LP::PivType::Frac;

            try {
                TG = Graph::TourGraph(tour_edges, edges, perm);
            } CMR_CATCH_PRINT_THROW("updating tour graph", err);
            
            continue;
        }

        Data::SupportGroup &supp_data = core_lp.supp_data;
        Sep::Separator separator(graph_data, best_data, supp_data,
                                  karp_part);

        try {
            if (!separator.find_cuts(TG)) {
                report_piv(piv, round);
                cout << "\tNo cuts found.\n";
                break;
            }
        } CMR_CATCH_PRINT_THROW("finding cuts", err);

        try {
            core_lp.pivot_back();
        } CMR_CATCH_PRINT_THROW("pivoting back", err);

        try {
            core_lp.add_cuts(separator.segment_q());
            core_lp.add_cuts(separator.fastblossom_q());
            core_lp.add_cuts(separator.blockcomb_q());
            core_lp.add_cuts(separator.connect_cuts_q());

            core_lp.add_cuts(separator.simpleDP_q());
        } CMR_CATCH_PRINT_THROW("adding cuts", err);
    }

    timer.stop();
    timer.report(false);
    cout << "\tFinal LP has " << core_lp.num_rows() << " rows, "
         << core_lp.num_cols() << " cols.\n";
    
    int stcount = 0;
    int dpcount = 0;
    
    for (const Sep::HyperGraph &H : core_lp.ext_cuts.get_cuts())
        if (H.cut_type() != CutType::Domino)
            ++stcount;
        else
            ++dpcount;
    cout << "\t" << stcount << " standard cuts, " << dpcount << " dp cuts.\n";
    cout << "\n";
    return piv;
}

}
