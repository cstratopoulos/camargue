#ifndef CMR_SOLVER_H
#define CMR_SOLVER_H

#include "core_lp.hpp"
#include "karp.hpp"
#include "datagroups.hpp"
#include "pricer.hpp"
#include "brancher.hpp"
#include "util.hpp"

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
    LP::PivType cutting_loop(bool do_price, bool try_recover);

    /// Embed cutting_loop in an augment and branch and cut search.
    LP::PivType abc(bool do_price);

    const Data::Instance &inst_info() const
        { return tsp_instance; } /// Get the Instance being used.
    
    const Data::GraphGroup &graph_info() const
        { return graph_data; } /// Get the GraphGroup being used.

    const Data::BestGroup &best_info() const
        { return best_data; } /// Get the BestGroup of the best tour found.

    const LP::TourBasis &tour_basis() const
        { return core_lp.tour_base; } /// Get the basis of the best tour found.

    const LP::CoreLP &get_core_lp() const
        { return core_lp; } /// Get the active CoreLP relaxation.

    /// Which separation routines should be called.
    struct CutSel {
        bool segment = true; //!< Primal SECs.
        bool fast2m = true; //!< Standard fast blossom heuristics.
        bool blkcomb = true; //!< Standard block comb heuristics.
        bool simpleDP = true; //!< Primal simple DP separation.
        bool safeGMI = false; //!< Primal safe Gomory cuts.
        bool connect = true; //!< Standard connect cut SECs.
    };

    CutSel cut_sel;

private:
    void report_piv(CMR::LP::PivType piv, int round, int num_pruned,
                    bool full_opt);

    LP::PivType cut_and_piv(int &round, int &num_pruned, bool do_price);

    LP::PivType frac_recover();
    
    Data::Instance tsp_instance;
    Data::KarpPartition karp_part;
    Data::GraphGroup graph_data;
    Data::BestGroup best_data;

    LP::CoreLP core_lp;

    Graph::TourGraph TG;

    std::unique_ptr<Price::Pricer> edge_pricer;
    std::unique_ptr<ABC::Brancher> brancher;

    OutPrefs output_prefs;
};

}

#endif
