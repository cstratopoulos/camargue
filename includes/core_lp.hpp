/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief CORE LP RELAXATIONS OF TSP INSTANCES
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CORE_LP_H
#define CMR_CORE_LP_H

#include "lp_interface.hpp"
#include "datagroups.hpp"

#include <vector>

namespace CMR {
namespace LP {

/** Class for storing the core lp associated to a TSP instance and pivoting.
 * This class contains the edges and constraints currently under consideration
 * in an lp relaxation of a TSP instance. 
 */
class CoreLP : public Relaxation {
public:
    CoreLP(CMR::Data::GraphGroup &graph_data_,
           CMR::Data::BestGroup &best_data_);

    CMR::LP::PivType primal_pivot();
    void pivot_back();

private:
    void handle_aug();
    void handle_fathom();
    
    CMR::Data::GraphGroup &graph_data;
    CMR::Data::BestGroup &best_data;
    CMR::Data::SupportGroup supp_data;

    struct TourBasis {
        TourBasis() = default;
        TourBasis(CMR::Graph &graph, CMR::Data::BestGroup &best_data);

        std::vector<double> best_tour_edges;
     
        std::vector<int> colstat;
        std::vector<int> rowstat;
    } tour_base;
    
    std::vector<double> lp_edges;
    std::vector<double> feas_stat;
};

}
}

#endif
