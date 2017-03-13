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
#include "dfs_brancher.hpp"

#if CMR_HAVE_SAFEGMI

#include "safeGMI.hpp"

#endif

#include "pricer.hpp"
#include "util.hpp"

#include <functional>
#include <memory>
#include <string>

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
        bool safeGMI = false; //!< Primal safe Gomory cuts.
        bool connect = true; //!< Standard connect cut SECs.
    } cut_sel;

private:
    /// Print info about \p piv and the size of the CoreLP.
    void report_lp(LP::PivType piv);

    void report_cuts(); //!< Report the number and types of cuts in CoreLP.

    void report_aug(bool piv_aug); //!< Output info about a new tour found.

    void initial_prints(); //!< Handles writing initial data to file.

    /// Should cut_and_piv start again with the easier separation routines.
    bool restart_loop(LP::PivType piv, double delta_metric);

    /// Should cut_and_piv return \p piv for augmentation or pricing.
    bool return_pivot(LP::PivType piv);

    /// A loop of primal pivoting and cutting plane generation.
    LP::PivType cut_and_piv(bool do_price);

    LP::PivType abc_dfs(bool do_price);

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

    std::unique_ptr<ABC::DFSbrancher> dfs_brancher;

    bool branch_engaged = false; //!< Is an ABC search active.

    OutPrefs output_prefs;

    int num_augs = 0;
};

}

#endif
