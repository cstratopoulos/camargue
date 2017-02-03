#include "config.hpp"
#include "solver.hpp"
#include "separator.hpp"

#if CMR_HAVE_SAFEGMI

#include "safeGMI.hpp"

#endif

#include "timer.hpp"
#include "err_util.hpp"

#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>

using std::abs;
using std::ceil;

using std::cout;
using std::cerr;
using std::endl;

using std::string;

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::unique_ptr;
using std::vector;
using std::function;


namespace CMR {

using CutType = Sep::HyperGraph::Type;
using PivType = LP::PivType;
namespace Eps = Epsilon;

template<class Qtype>
bool call_separator(const function<bool()> &sepcall, const Qtype &sep_q,
                    PivType &piv, LP::CoreLP &core_lp,
                    const double tourlen, double &prev_val,
                    double &total_delta,double &delta_ratio, int &num_pruned)
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

void Solver::report_piv(PivType piv, int round, int num_pruned, bool full_opt)
{
    int rowcount = core_lp.num_rows();
    cout << "\n\t\tRound " << round << "\n"
         << "\tLP objective value: " << core_lp.get_objval()
         << ", dual feasible: " << core_lp.dual_feas() << "\n";
    cout << "\t" << rowcount << " rows, "
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
        cout << "\tAugmented to new tour, pruned "
             << num_pruned << " slack cuts from LP\n";
        break;
    default:
        cout << "\t Pivot status: \t" << piv << "\n";
    }
    
    cout << endl;
}

PivType Solver::cutting_loop(bool do_price, bool try_recover)
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


    CMR::Timer timer(tsp_instance.problem_name());
    
    timer.start();    

    while (true) {
        ++auground;
        int num_pruned = 0;

        try {
            piv = cut_and_piv(round, num_pruned, do_price);
        } CMR_CATCH_PRINT_THROW("invoking cut and piv", err);

        if (piv == PivType::FathomedTour) {
            report_piv(piv, round, num_pruned, !do_price);
            
            if (do_price) {
                try {
                    if (edge_pricer->gen_edges(piv) == Price::ScanStat::Full) {
                        core_lp.rebuild_basis();
                        continue;
                    } else
                        break;
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }            
            break;
        }

        if (piv == PivType::Tour) {
            report_piv(piv, round, num_pruned, false);

            if (do_price) {
                try {
                    edge_pricer->gen_edges(piv);
                    core_lp.rebuild_basis();
                } CMR_CATCH_PRINT_THROW("adding edges to core", err);
            }
            
            try {
                TG = Graph::TourGraph(tour_edges, edges, perm);
            } CMR_CATCH_PRINT_THROW("updating tour graph", err);
            
            continue;
        }

        if (try_recover)
            try {
                int prev_rows = core_lp.num_rows();
                piv = frac_recover();
                num_pruned = prev_rows - core_lp.num_rows();
                if (piv == PivType::Tour) {
                    report_piv(piv, round, num_pruned, false);
                    TG = Graph::TourGraph(tour_edges, edges, perm);

                    continue;
                }
            } CMR_CATCH_PRINT_THROW("trying to recover from frac tour", err);

        report_piv(piv, round,  0, false);
        cout << "\tNo cuts found.\n";
        break;
    }

    timer.stop();
    timer.report(true);

    if (piv == PivType::FathomedTour && do_price)
        report_piv(piv, round, 0, true);

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

    cout << "\t" << subcount << " SECs, " << combcount << " combs/blossoms, "
         << dpcount << " dp cuts, " << gmicount << " GMI cuts.\n";
    cout << "\n";
    return piv;
}

PivType Solver::abc(bool do_price)
{
    runtime_error err("Problem in Solver::abc");

    PivType piv = cutting_loop(do_price, true);
    
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
        
        try {
            core_lp.factor_basis();
            piv = cutting_loop(do_price, false);
        } CMR_CATCH_PRINT_THROW("cutting on branch prob", err);
        
        cout << "\tTASK STATUS: ";
        if (piv == PivType::FathomedTour) {
            brancher->pop_problem(ProbStat::Pruned);
        } else if (piv == PivType::Frac) {
            brancher->pop_problem(ProbStat::Seen);
        } else {
            cerr << "Pivot status " << piv << " in abc loop.\n";
            throw runtime_error("Invalid piv stat for further branching.");
        }
        
        cout << "\tCalling next prob....\n";
        prob = brancher->next_prob();
    }

    cout << "\n\t\tABC search complete.\n\n";
    
    return piv;    
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
                else
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
