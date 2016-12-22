#ifndef CMR_SOLVER_H
#define CMR_SOLVER_H

#include "core_lp.hpp"
#include "datagroups.hpp"
#include "separator.hpp"
#include "util.hpp"

#include <string>

namespace CMR {

class Solver {
public:
    /** Construct a Solver from a TSPLIB instance. */
    Solver(const std::string &tsp_fname, const int seed,
           const CMR::OutPrefs outprefs);

    /** Construct Solver from TSPLIB instance with starting tour from file. */
    Solver(const std::string &tsp_fname, const std::string &tour_fname,
           const int seed, const CMR::OutPrefs outprefs);

    /** Construct Solver from randomly generated Euclidean TSP instance. */
    Solver(const int seed, const int node_count, const int gridsize,
           const CMR::OutPrefs outprefs);

    CMR::LP::PivType cutting_loop();

    CMR::LP::Relaxation &relax() { return core_lp; }

    const CMR::LP::TourBasis &tour_basis() const { return core_lp.tour_base; }

    const Data::BestGroup &best_info() const { return best_data; }

    

private:
    void report_piv(CMR::LP::PivType piv, int round);
    
    Data::Instance tsp_instance;
    Data::KarpPartition karp_part;
    Data::GraphGroup graph_data;
    Data::BestGroup best_data;

    CMR::LP::CoreLP core_lp;

    CMR::OutPrefs output_prefs;
};

}

#endif
