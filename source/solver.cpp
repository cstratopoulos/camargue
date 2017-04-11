/**
 * @file
 * @brief Implementation of higher level or public Solver methods.
 */

#include "solver.hpp"
#include "io_util.hpp"
#include "timer.hpp"
#include "err_util.hpp"
#include "config.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

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

Solver::Solver(const string &fname, int seed, Graph::EdgePlan eplan,
               OutPrefs outprefs)
try : tsp_instance(fname, make_seed(seed)),
      karp_part(tsp_instance),
      core_graph(tsp_instance, eplan), best_data(tsp_instance, core_graph),
      core_lp(core_graph, best_data),
      output_prefs(outprefs)
{
    initial_prints();
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Solver TSPLIB constructor failed.");
}

Solver::Solver(const std::string &tsp_fname, int seed, OutPrefs outprefs)
    : Solver(tsp_fname, seed, Graph::EdgePlan::Linkern, outprefs) {}

Solver::Solver(const string &fname, const string &tour_fname,
               int seed, Graph::EdgePlan eplan, OutPrefs outprefs)
try : tsp_instance(fname, make_seed(seed)),
      karp_part(tsp_instance),
      core_graph(tsp_instance, eplan),
      best_data(tsp_instance, core_graph, tour_fname),
      core_lp(core_graph, best_data),
      output_prefs(outprefs)
{
    initial_prints();
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Solver TSPLIB/tour constructor failed.");
}

Solver::Solver(const std::string &tsp_fname, const std::string &tour_fname,
               int seed, OutPrefs outprefs)
    : Solver(tsp_fname, tour_fname, seed, Graph::EdgePlan::Linkern, outprefs)
{}

Solver::Solver(int seed, int node_count, int gridsize, Graph::EdgePlan eplan,
               OutPrefs outprefs)
try : tsp_instance(make_seed(seed), node_count, gridsize),
      karp_part(tsp_instance),
      core_graph(tsp_instance, eplan), best_data(tsp_instance, core_graph),
      core_lp(core_graph, best_data),
      output_prefs(outprefs)
{
    initial_prints();
} catch (const exception &e) {
    cerr << e.what() << endl;
    throw runtime_error("Solver random constructor failed.");
}

Solver::Solver(int seed, int node_count, int gridsize, OutPrefs outprefs)
    : Solver(seed, node_count, gridsize, Graph::EdgePlan::Linkern, outprefs) {}

void Solver::choose_cuts(CutSel::Presets preset)
{
    using CutPre = CutSel::Presets;
    cut_sel = CutSel{};
    if (preset == CutPre::Vanilla)
        return;

    if (preset == CutPre::Aggressive || preset == CutPre::Sparse) {
        cut_sel.localcuts = true;
        cut_sel.decker = true;
        cut_sel.handling = true;
        cut_sel.teething = true;
        cut_sel.consec1 = true;
        cut_sel.tighten = true;
        cut_sel.tighten_pool = true;
        if (preset == CutPre::Sparse)
            cut_sel.safeGMI = true;
        return;
    }

    throw logic_error("Unimplemented CutSel::Presets choice");
}

void Solver::report_lp(PivType piv)
{
    int rowcount = core_lp.num_rows();
    cout << "\tPivot status:\t" << piv << endl
         << "\tLP objective value: " << core_lp.get_objval()
         << ", dual feasible: " << core_lp.dual_feas() << endl;
    cout << "\t" << rowcount << " rows, "
         << core_lp.num_cols() << " cols in LP.\n" << endl;
}

ostream &operator<<(ostream &os, Solver::Aug aug)
{
    switch (aug) {
    case Solver::Aug::Init:
        os << "Initial Tour";
        break;
    case Solver::Aug::Piv:
        os << "Primal Pivot";
        break;
    case Solver::Aug::Xtour:
        os << "X-tour Heuristic";
        break;
    case Solver::Aug::Branch:
        os << "Branch Tour";
        break;
    }

    return os;
}

/**
 * Prints the progress in tour improvement.
 * @param aug_type the Solver::Aug describing the mode of augmentation.
 */
void Solver::report_aug(Aug aug_type)
{
    if (output_prefs.prog_bar)
        cout << "\n";
    cout << "\tTour " << ++num_augs << ": "
         << static_cast<int>(best_data.min_tour_value) << ", "
         << aug_type << endl;

    aug_chart.emplace_back(aug_type, core_lp.get_objval());

    if (!output_prefs.save_tour_edges && !output_prefs.save_tour)
        return;

    string infix = file_infix();

    if (output_prefs.save_tour || output_prefs.save_tour_edges){
        cout << "\tWrote tour ";
        if (output_prefs.save_tour) {
            string tour_fname = tsp_instance.problem_name() + infix + "sol";
            util::write_tour_nodes(best_data.best_tour_nodes, tour_fname);
            cout << "nodes to " << tour_fname;
            if (output_prefs.save_tour_edges)
                cout << ", ";
        }

        if (output_prefs.save_tour_edges) {
            string edges_fname = tsp_instance.problem_name() + "_tour" + infix +
            "x";
            util::write_tour_edges(best_data.best_tour_edges,
                                   core_graph.get_edges(),
                                   core_graph.node_count(), edges_fname);
            cout << "edges to " << edges_fname;
        }
        cout << endl;
    }
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
         << " combs/blossoms/clique things, " << dpcount << " dp cuts, "
         << gmicount << " GMI cuts. \n\t("
         << (core_lp.num_rows() - tsp_instance.node_count())
         << " cuts total, " << core_lp.num_cols() << " cols, "
         << core_lp.ext_cuts.pool_count() << " cuts in pool)\n" << endl;
}

void Solver::initial_prints()
{
    time_overall.start();
    time_piv.start(); time_piv.stop();
    time_price.start(); time_price.stop();
    time_branch.start(); time_branch.stop();

    for (auto &kv : sep_times) {
        kv.second.first.start();
        kv.second.first.stop();
    }

    bool want_xy = output_prefs.dump_xy;
    bool want_tour = output_prefs.save_tour;
    bool want_edges = output_prefs.save_tour_edges;

    string pname = tsp_instance.problem_name();

    aug_chart.emplace_back(Aug::Init, best_data.min_tour_value);

    if (want_xy) {
        if (tsp_instance.ptr()->x == NULL)
            cerr << "Dumping XY coords incompatible with " << pname << endl;
        else {
            string xy_fname = pname + ".xy";
            util::write_xy_coords(tsp_instance.ptr()->x,
                                  tsp_instance.ptr()->y,
                                  tsp_instance.node_count(),
                                  xy_fname);
            cout << "Wrote XY coords to " << xy_fname << endl;
        }
    }

    string infix = file_infix();

    if (want_tour) {
        string tour_fname = pname + infix + "sol";
        util::write_tour_nodes(best_data.best_tour_nodes,
                               tour_fname);
        cout << "Wrote starting tour to " << tour_fname << endl;
    }

    if (want_edges) {
        string edges_fname = pname + "_tour" + infix + "x";
        util::write_tour_edges(best_data.best_tour_edges,
                               core_graph.get_edges(),
                               tsp_instance.node_count(), edges_fname);
        cout << "Wrote starting tour edges to " << edges_fname << endl;
    }
}

std::string Solver::file_infix()
{
    if (!output_prefs.gif_tour) {
        return ".";
    } else {
        return "." + std::to_string((int) best_data.min_tour_value) + ".";
    }
}

void Solver::set_lowerbound(double lb)
{
    target_lb = std::ceil(lb);
    cout << "Set target lower bound to " << target_lb << endl;
}

bool Solver::lb_fathom()
{
    return best_data.min_tour_value <= target_lb;
}


PivType Solver::cutting_loop(bool do_price, bool try_recover, bool pure_cut)
{
    runtime_error err("Problem in Solver::cutting_loop");

    if (lb_fathom()) {
        cout << "Starting tour already matches target LB, returning optimal."
             << endl;
        return PivType::FathomedTour;
    }

    if (do_price)
        try {
            time_price.resume();
            edge_pricer = util::make_unique<Price::Pricer>(core_lp,
                                                           tsp_instance,
                                                           core_graph);
            time_price.stop();
            edge_pricer->verbose = output_prefs.verbose;
        } CMR_CATCH_PRINT_THROW("instantiating/allocating Pricer", err);

    core_lp.verbose = output_prefs.verbose;

    PivType piv = PivType::Frac;
    bool elim_during = true;

    int rounds = 0;

    while (!lb_fathom()) {
        ++rounds;
        if (output_prefs.verbose)
            cout << "|| Cutting loop pass " << rounds << endl;

        try {
            piv = cut_and_piv(do_price);
        } CMR_CATCH_PRINT_THROW("invoking cut and piv", err);

        if (piv == PivType::Subtour && cut_sel.connect)
            throw logic_error("Left cut and piv with integral subtour");

        if (piv == PivType::FathomedTour) {
            if (do_price) {
                if (pure_cut)
                    cout << "\tTour optimal for edge set...\n";
                edge_pricer->verbose = output_prefs.verbose + pure_cut;
                try {
                    time_price.resume();
                    Price::ScanStat scan_stat =
                    edge_pricer->gen_edges(piv, elim_during &&pure_cut);
                    time_price.stop();

                    if (scan_stat == Price::ScanStat::Full) {
                        time_piv.resume();
                        core_lp.pivot_back(false);
                        time_piv.stop();
                        continue;
                    } else
                        break;
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }
            break;
        }

        if (piv == PivType::Tour) {
            if (core_lp.active_tourlen() < best_data.min_tour_value) {
                try { core_lp.active_tour.best_update(best_data); }
                CMR_CATCH_PRINT_THROW("pivot updating best data", err);

                report_aug(Aug::Piv);
                if (lb_fathom()) {
                    cout << "Augmenting pivot matches target LB, "
                         << "returning optimal." << endl;
                    piv = PivType::FathomedTour;
                    continue;
                }
            }

            if (do_price) {
                try {
                    time_price.resume();
                    edge_pricer->gen_edges(piv, false);
                    time_price.stop();
                }
                CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }

            continue;
        }


        if (try_recover) {
            try {
                if (frac_recover() == PivType::Tour) {
                    piv = PivType::Tour;
                    if (core_lp.active_tourlen() < best_data.min_tour_value) {
                        core_lp.active_tour.best_update(best_data);

                        report_aug(Aug::Xtour);
                        if (lb_fathom()) {
                            cout << "X-tour pivot matches target LB, "
                                 << "returning optimal." << endl;
                            piv = PivType::FathomedTour;
                            continue;
                        }
                    }
                    continue;
                }
            } CMR_CATCH_PRINT_THROW("trying to recover from frac tour", err);
        }

        if (pure_cut) {
            cout << "\tNo cuts found.\n";
        }
        break;
    }

    if (pure_cut) {
        cout << "\tPivot status " << piv << ", obj val "
             << core_lp.get_objval() << endl;
        report_cuts();
    }
    return piv;
}

Solver::~Solver()
{
    if (!output_prefs.detailed_stats)
        return;
    bool want_cpu = false;
#ifdef CMR_USE_OMP
    want_cpu = true;
#endif

    time_overall.stop();
    cout << "\n";
    time_piv.report(false);
    time_price.report(false);
    if (branch_engaged)
        time_branch.report(false);
    cout << "\n";
    for(const auto &kv : sep_times) {
        if (kv.second.second == false)
            continue;

        const string &name = kv.first;
        const Timer &T = kv.second.first;

        if (name == "ExactBlossoms" || name == "SimpleDP")
            T.report(want_cpu);
        else
            T.report(false);
    }

    cout << "\n";

    time_overall.report(want_cpu);
}

}
