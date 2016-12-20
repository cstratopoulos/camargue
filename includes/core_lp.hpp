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

#include <vector>

namespace CMR {

class Solver;

namespace LP {

struct TourBasis {
    TourBasis() = default;
    TourBasis(CMR::Graph &graph, CMR::Data::BestGroup &best_data);

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
    CoreLP(CMR::Data::GraphGroup &graph_data_,
           CMR::Data::BestGroup &best_data_);

    double opt();

    CMR::LP::PivType primal_pivot();
    void pivot_back();

    void add_cuts(CMR::Sep::LPcutList &cutq);
    void add_cuts(CMR::Sep::CutQueue<CMR::Sep::dominoparity> &dp_q);


    CMR::Data::SupportGroup supp_data;

    friend class CMR::Solver;

private:
    void handle_aug();
    void handle_fathom();
    
    CMR::Data::GraphGroup &graph_data;
    CMR::Data::BestGroup &best_data;

    TourBasis tour_base;
    
    std::vector<double> lp_edges;
    std::vector<double> feas_stat;
};

}
}

#endif
