/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief TSP Solver class header.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_SOLVER_H
#define CMR_SOLVER_H

#include "core_lp.hpp"
#include "karp.hpp"
#include "datagroups.hpp"
#include "separator.hpp"
#include "meta_sep.hpp"
#include "abc_nodesel.hpp"

#include "pricer.hpp"
#include "err_util.hpp"
#include "util.hpp"
#include "timer.hpp"

#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <stdexcept>

namespace CMR {

namespace Sep {
class SafeGomory;
}

/// Solution of TSP instances.
class Solver {
public:
    /**@name TSPLIB constructors, each with or without specifying EdgePlan. */
    ///@{
    Solver(const std::string &tsp_fname, int seed, Graph::EdgePlan eplan,
           OutPrefs outprefs);
    Solver(const std::string &tsp_fname, int seed, OutPrefs outprefs);

    ///Specify a starting tour from file.
    Solver(const std::string &tsp_fname, const std::string &tour_fname,
           int seed, Graph::EdgePlan eplan, OutPrefs outprefs);
    Solver(const std::string &tsp_fname, const std::string &tour_fname,
           int seed, OutPrefs outprefs);
    ///@}

    /**@name Random Euclidean constructors, with or without EdgePlan. */
    ///@{
    Solver(int seed, int node_count, int gridsize, Graph::EdgePlan eplan,
           OutPrefs outprefs);
    Solver(int seed, int node_count, int gridsize, OutPrefs outprefs);
    ///@}

    ~Solver();

    void set_lowerbound(double lb); //!< Set a target for early termination.

    /// Run a primal cutting plane loop of pivoting and cut generation.
    LP::PivType cutting_loop(bool do_price, bool try_recover,
                             bool pure_cut);

    template <typename SelectionRule>
    LP::PivType abc(bool do_price); //!< Augment-branch-cut search.

    const Data::Instance &inst_info() const { return tsp_instance; }
    const Graph::CoreGraph &graph_info() const { return core_graph; }
    const Data::BestGroup &best_info() const{ return best_data; }
    const LP::ActiveTour &active_tour() const { return core_lp.active_tour; }
    const LP::CoreLP &get_core_lp() const { return core_lp; }

    /// Which separation routines should be called.
    struct CutSel {
        /// Preset selection routine choices. Each contains the one before it.
        enum Presets {
            Vanilla, //!< The default initializations indicated below.
            Aggressive, //!< Add cut metamorphoses/local cuts.
            Sparse, //!< Add safe Gomory cuts.
        };

        bool cutpool = true; //!< Cuts from a cut pool.

        /**@name Primal template cuts. */
        ///@{
        bool segment = true; //!< Exact SECs.
        bool ex2m = true; //!< Exact blossoms.
        bool simpleDP = true; //!< Partitioned simple DP cuts.
        ///@}

        /**@name Fast standard heuristics. */
        ///@{
        bool fast2m = true; //!< Odd component (PH/GH) blossoms.
        bool blkcomb = true; //!< Block combs.
        bool connect = true; //!< Connected component SECs.
        ///@}

        /**@name Standard non-template cuts and metamorphoses. */
        ///@{
        bool localcuts = false;
        bool decker = false;
        bool handling = false;
        bool teething = false;

        bool consec1 = false;

        bool tighten = false;
        bool tighten_pool = false;
        ///@}

        bool safeGMI = false; //!< Primal safe Gomory cuts.
    } cut_sel;

    void choose_cuts(CutSel::Presets preset); //!< Set the choice of cuts.

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
    void report_lp(LP::PivType piv); //!< Report on \p piv and core_lp.
    void report_cuts(); //!< Report the number and types of cuts in CoreLP.
    void report_aug(Aug aug_type); //!< Output info about a new tour found.
    void initial_prints(); //!< Handles writing initial data to file.

    std::string file_infix(); //!< Returns an infix for data saved to file.

    /// Compute the Padberg-Hong-esque delta ratio.
    static double PH_delta(double new_val, double prev_val, double tourlen);

    /// Should cut_and_piv start again with the easier separation routines.
    bool restart_loop(LP::PivType piv, double delta_metric);

    /// Should cut_and_piv return \p piv for augmentation or pricing.
    bool return_pivot(LP::PivType piv);

    bool lb_fathom(); //!< Does the current best tour match or beat target_lb.

    /// A loop of primal pivoting and cutting plane generation.
    LP::PivType cut_and_piv(bool do_price);

    /// Check if \p prob is prunable based on an optimize/price estimate.
    void opt_check_prunable(bool do_price, ABC::BranchNode &prob);

    LP::PivType abc_bcp(bool do_price); //!< The branch-cut-price loop.

    LP::PivType frac_recover(); //!< Use x-tour heuristic on fractional pivot.

    /**@name Separator resetting routines.
     * These methods reset their argument with data about the current pivot.
     */
    ///@{
    void reset_separator(std::unique_ptr<Sep::Separator> &S);
    void reset_separator(std::unique_ptr<Sep::MetaCuts> &MS);

    void reset_separator(std::unique_ptr<Sep::SafeGomory> &GS);
    ///@}

    /// Tracking objective values of pivots within a cut_and_piv loop.
    struct PivStats {
        /// Construct PivStats from the first pivot before adding cuts.
        PivStats(double first_piv) :
            initial_piv(first_piv), prev_val(first_piv), lowest_piv(first_piv)
            {}

        void update(double new_val, double tourlen);

        void report_extrema();

        const double initial_piv; //!< The first pivot val before adding cuts.
        double prev_val = 0.0; //!< The previous pivot val.
        double delta_ratio = 0.0; //!< The Padberg-Hong delta ratio.
        double lowest_piv = 0.0; //!< The lowest piv val seen.
        double max_ratio = 0.0; //!< The highest delta ratio seen.
        double first_last_ratio = 0.0; //!< Delta ratio for first and last piv.

        bool found_cuts = false; //!< Have any cuts been found.
    };

    /// Method template for calling a separation routine in cut_and_piv.
    template <typename Qtype>
    bool call_separator(const std::function<bool()> &sepcall,
                        Qtype &sep_q, const std::string sep_name,
                        LP::PivType &piv, PivStats &piv_stats,
                        bool pivback_prune);

    double target_lb{-std::numeric_limits<double>::max() + 1.0};

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

    std::array<char, 80> p_bar;
    void place_pivot(double low_limit, double best_tourlen, double piv_val);

    Timer time_overall{"Solver Overall"};
    Timer time_piv{"Pivoting", &time_overall};
    Timer time_price{"Pricing", &time_overall};
    Timer time_branch{"Branching", &time_overall};

    using TimerCalled = std::pair<Timer, bool>;

    std::map<std::string, TimerCalled> sep_times{
        {"CutPool", {Timer("CutPool", &time_overall), false}},
        {"SegmentCuts", {Timer("SegmentCuts", &time_overall), false}},
        {"ConnectCuts", {Timer("ConnectCuts", &time_overall), false}},
        {"ExactSub", {Timer("ExactSub", &time_overall), false}},
        {"FastBlossoms", {Timer("FastBlossoms", &time_overall), false}},
        {"ExactBlossoms", {Timer("FastBlossoms", &time_overall), false}},
        {"BlockCombs", {Timer("BlockCombs", &time_overall), false}},
        {"SimpleDP", {Timer("SimpleDP", &time_overall), false}},
        {"LocalCuts", {Timer("LocalCuts", &time_overall), false}},
        {"Decker", {Timer("Decker", &time_overall), false}},
        {"Handling", {Timer("Handling", &time_overall), false}},
        {"Teething", {Timer("Teething", &time_overall), false}},
        {"Consec1", {Timer("Consec1", &time_overall), false}},
        {"Tighten", {Timer("Tighten", &time_overall), false}},
        {"TightenPool", {Timer("TightenPool", &time_overall), false}},
        {"SafeGMI", {Timer("SafeGMI", &time_overall), false}},
    };
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
            time_price.resume();
            edge_pricer->elim_edges(true);
            core_lp.primal_opt();
            time_price.stop();
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

    try {
        branch_controller = util::make_unique<SelectionRule>(tsp_instance,
                                                             best_info(),
                                                             graph_info(),
                                                             core_lp);
    } CMR_CATCH_PRINT_THROW("allocating/instantiating Brancher", err);

    branch_engaged = true;
    branch_controller->verbose = output_prefs.verbose;

    if (cut_sel.safeGMI) {
        cout << "(Disabling GMI and purging cuts for branching.....)\n";
        cut_sel.safeGMI = false;
        try { core_lp.purge_gmi(true); }
        CMR_CATCH_PRINT_THROW("dumping gmi cuts before abc", err);
    }

    try { piv = abc_bcp(do_price); }
    CMR_CATCH_PRINT_THROW("running abc_bcp", err);

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

    return piv;
}


}

#endif
