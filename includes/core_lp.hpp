/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief CORE LP RELAXATIONS OF TSP INSTANCES
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CORE_LP_H
#define CMR_CORE_LP_H

#include "lp_interface.hpp"
#include "cc_lpcuts.hpp"
#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "util.hpp"

#include <vector>

namespace CMR {

class Solver;

namespace LP {

struct TourBasis {
    TourBasis() = default;
    TourBasis(const Graph::CoreGraph &graph,
              const Data::BestGroup &best_data);

    std::vector<double> best_tour_edges;
     
    std::vector<int> colstat;
    std::vector<int> rowstat;
};

/** Class for storing the core lp associated to a TSP instance and pivoting.
 * This class contains the edges and constraints currently under consideration
 * in an lp relaxation of a TSP instance. 
 */
class CoreLP : public Relaxation {
public:
    CoreLP(Data::GraphGroup &graph_data_,
           Data::BestGroup &best_data_);

    LP::PivType primal_pivot();
    void pivot_back();

    void add_cuts(const Sep::LPcutList &cutq);
    void add_cuts(const Sep::CutQueue<Sep::dominoparity> &dp_q);

    void add_edges(const std::vector<Graph::Edge> &add_batch);

    const Sep::ExternalCuts &external_cuts() const
        { return ext_cuts; } //<! Const reference to the external cut reps.

    friend class CMR::Solver;

private:
    void handle_aug();

    void rebuild_basis();
    
    Data::GraphGroup &graph_data;
    Data::BestGroup &best_data;
    Data::SupportGroup supp_data;

    Sep::ExternalCuts ext_cuts;

    TourBasis tour_base;
    
    std::vector<double> lp_edges;
    std::vector<double> feas_stat;


};

}
}

#endif
