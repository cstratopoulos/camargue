#ifndef CMR_SOLVER_H
#define CMR_SOLVER_H

#include "datagroups.hpp"
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

private:
    CMR::Data::Instance tsp_instance;
    CMR::Data::GraphGroup graph_data;
    CMR::Data::BestGroup best_data;

    CMR::OutPrefs output_prefs;
};

}

#endif
