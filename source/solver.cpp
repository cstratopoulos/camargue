#include "solver.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

using std::string;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

inline static int make_seed(const int seed)
{
    return (seed > 0) ? seed : (int) real_zeit();
}

Solver::Solver(const string &fname, const int seed, const OutPrefs outprefs)
try : tsp_instance(fname, make_seed(seed)),
      graph_data(tsp_instance), best_data(tsp_instance, graph_data),
      output_prefs(outprefs)
    
{
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver TSPLIB constructor failed.");
}

Solver::Solver(const string &fname, const string &tour_fname,
               const int seed, const OutPrefs outprefs)
try : tsp_instance(fname, make_seed(seed)),
      graph_data(tsp_instance),
      best_data(tsp_instance, graph_data, tour_fname),
      output_prefs(outprefs)
{
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver TSPLIB/tour constructor failed.");
}

Solver::Solver(const int seed, const int node_count,
               const int gridsize, const OutPrefs outprefs)
try : tsp_instance(make_seed(seed), node_count, gridsize),
      graph_data(tsp_instance), best_data(tsp_instance, graph_data),
      output_prefs(outprefs)
{
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Solver random constructor failed.");
}

}
