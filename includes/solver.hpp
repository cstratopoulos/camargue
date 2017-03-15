/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief TSP Solver class header.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_SOLVER_H
#define CMR_SOLVER_H

#include "config.hpp"
#include "core_lp.hpp"
#include "karp.hpp"
#include "datagroups.hpp"
#include "separator.hpp"
#include "pool_sep.hpp"
#include "abc_nodesel.hpp"

#if CMR_HAVE_SAFEGMI

#include "safeGMI.hpp"

#endif

#include "pricer.hpp"
#include "err_util.hpp"
#include "util.hpp"
#include "timer.hpp"

#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>

namespace CMR {

/// Solution of TSP instances.
class Solver {
public:
    /// Construct a Solver from a TSPLIB instance.
    Solver(const std::string &tsp_fname, const int seed,
           const OutPrefs outprefs);

    /// Construct Solver from TSPLIB instance with starting tour from file.
    Solver(const std::string &tsp_fname, const std::string &tour_fname,
           const int seed, const OutPrefs outprefs);

    /// Construct Solver from randomly generated Euclidean TSP instance.
    Solver(const int seed, const int node_count, const int gridsize,
           const OutPrefs outprefs);

    /// Run a primal cutting plane loop of pivoting and cut generation.
    LP::PivType cutting_loop(bool do_price, bool try_recover,
                             bool pure_cut);

    /// Embed cutting_loop in an augment and branch and cut search.
    template <typename SelectionRule>
    LP::PivType abc(bool do_price);

    const Data::Instance &inst_info() const { return tsp_instance; }

    const Graph::CoreGraph &graph_info() const { return core_graph; }

    const Data::BestGroup &best_info() const{ return best_data; }

    const LP::ActiveTour &active_tour() const { return core_lp.active_tour; }

    const LP::CoreLP &get_core_lp() const { return core_lp; }

    /// Which separation routines should be called.
    struct CutSel {
        bool cutpool = true; //!< Cuts from a cut pool.
        bool segment = true; //!< Primal SECs.
        bool fast2m = true; //!< Standard fast blossom heuristics.
        bool blkcomb = true; //!< Standard block comb heuristics.
        bool ex2m = true; //!< Exact primal blossom separation.
        bool simpleDP = true; //!< Primal simple DP separation.
        bool localcuts = false; //!< Standard local cut separation.
        bool safeGMI = false; //!< Primal safe Gomory cuts.
        bool connect = true; //!< Standard connect cut SECs.
    } cut_sel;

    /// Types of augmentations that can take place.
    enum Aug : char {
        Init = 'I', //!< The starting tour.
        Piv = 'P', //!< An augmenting pivot.
        Xtour = 'X', //!< The x-tour/short LK heuristic.
        Branch = 'B', //!< An improving branch tour.
    };

    using AugObj = std::pair<Aug, double>;

    const std::vector<AugObj> &get_aug_chart() const { return aug_chart; }

private:
    /// Print info about \p piv and the size of the CoreLP.
    void report_lp(LP::PivType piv);

    void report_cuts(); //!< Report the number and types of cuts in CoreLP.

    void report_aug(Aug aug_type); //!< Output info about a new tour found.

    void initial_prints(); //!< Handles writing initial data to file.

    /// Should cut_and_piv start again with the easier separation routines.
    bool restart_loop(LP::PivType piv, double delta_metric);

    /// Should cut_and_piv return \p piv for augmentation or pricing.
    bool return_pivot(LP::PivType piv);

    /// A loop of primal pivoting and cutting plane generation.
    LP::PivType cut_and_piv(bool do_price);

    LP::PivType abc_bcp(bool do_price); //!< The branch-cut-price loop.

    /// Attempt to recover from a fractional pivot with the x-tour heuristic.
    LP::PivType frac_recover();

    /**@name Separator resetting routines.
     * These methods reset their argument with data about the current pivot.
     */
    ///@{

    void reset_separator(std::unique_ptr<Sep::Separator> &S);
    void reset_separator(std::unique_ptr<Sep::PoolCuts> &PS);

#if CMR_HAVE_SAFEGMI
    void reset_separator(std::unique_ptr<Sep::SafeGomory> &GS);
#endif

    ///@}

    /// Method template for calling a separation routine in cut_and_piv.
    template <typename Qtype>
    bool call_separator(const std::function<bool()> &sepcall,
                        Qtype &sep_q,
                        LP::PivType &piv, double &prev_val,
                        double &total_delta, double &delta_ratio,
                        double &lowest_piv);

    Data::Instance tsp_instance;
    Data::KarpPartition karp_part;
    Graph::CoreGraph core_graph;
    Data::BestGroup best_data;

    LP::CoreLP core_lp;

    std::unique_ptr<Price::Pricer> edge_pricer;

    std::unique_ptr<ABC::BaseBrancher> branch_controller;

    bool branch_engaged = false; //!< Is an ABC search active.

    OutPrefs output_prefs;

    int num_augs = 0;
    std::vector<AugObj> aug_chart;
};

std::ostream &operator<<(std::ostream &os, Solver::Aug aug);


///////////////////////// TEMPLATE IMPLEMENTATION /////////////////////////////


/**
 * @tparam SelectionRule the node selection to guide the ABC search. Should
 * be a derived class of ABC::BaseBrancher.
 * @param do_price will edge pricing be performed.
 */
template <typename SelectionRule>
LP::PivType Solver::abc(bool do_price)
{
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::runtime_error;
    using std::logic_error;

    using LP::PivType;

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

    if (do_price) {
        try {
            edge_pricer->elim_edges(true);
            core_lp.primal_opt();
            cout << "\tcol count " << core_lp.num_cols()
                 << ", opt objval " << core_lp.get_objval() << endl;
        } CMR_CATCH_PRINT_THROW("eliminating and optimizing", err);
    } else {
        try {
            core_lp.primal_opt();
        } CMR_CATCH_PRINT_THROW("optimizing at root", err);
        cout << "\tRoot LP optimized with obj val " << core_lp.get_objval()
             << endl;
    }


    cout << "\tCommencing ABC search....\n";
    cout << "Avg piv itcount " << core_lp.avg_itcount() << endl;
    Timer abct(tsp_instance.problem_name() + " ABC search");
    abct.start();

    try {
        branch_controller = util::make_unique<SelectionRule>(tsp_instance,
                                                             active_tour(),
                                                             best_info(),
                                                             graph_info(),
                                                             core_lp);
    } CMR_CATCH_PRINT_THROW("allocating/instantiating Brancher", err);

    branch_engaged = true;

    if (cut_sel.safeGMI) {
        cout << "(Disabling GMI and purging cuts for branching.....)\n";
        cut_sel.safeGMI = false;
        try { core_lp.purge_gmi(true); }
        CMR_CATCH_PRINT_THROW("dumping gmi cuts before abc", err);
    }

    try { piv = abc_bcp(do_price); }
    CMR_CATCH_PRINT_THROW("running abc_dfs", err);

    abct.stop();

    cout << "\n\tABC search completed, optimal tour has length "
         << best_data.min_tour_value << endl;


    int max_depth = 0;
    const ABC::BranchHistory &BH = branch_controller->get_history();
    for (const auto &B : BH)
        if (B.depth > max_depth)
            max_depth = B.depth;

    cout << "\t" << BH.size() << " branch nodes, max depth "
         << max_depth << endl;

    report_cuts();

    abct.report(true);

    return piv;
}


}

#endif
