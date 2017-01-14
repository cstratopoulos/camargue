#ifndef CMR_SOLVER_H
#define CMR_SOLVER_H

#include "core_lp.hpp"
#include "karp.hpp"
#include "datagroups.hpp"
#include "pricer.hpp"
#include "util.hpp"

#include <memory>
#include <string>

namespace CMR {

/// Class for solution of TSP instances.
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
    LP::PivType cutting_loop(bool do_price);

    /// Get the basis of the best tour found thus far.
    const LP::TourBasis &tour_basis() const { return core_lp.tour_base; }

    /// Get the active CoreLP relaxation.
    const LP::CoreLP &get_core_lp() const { return core_lp; }

    /// Get the BestGroup of the best tour found thus far.
    const Data::BestGroup &best_info() const { return best_data; }    

private:
    void report_piv(CMR::LP::PivType piv, int round);
    
    Data::Instance tsp_instance;
    Data::KarpPartition karp_part;
    Data::GraphGroup graph_data;
    Data::BestGroup best_data;

    LP::CoreLP core_lp;

    std::unique_ptr<Price::Pricer> edge_pricer;

    OutPrefs output_prefs;
};

}

#endif
