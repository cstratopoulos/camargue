/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Exact primal blossom separation. 
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BLOSSOMS_H
#define CMR_BLOSSOMS_H

#include "cut_structs.hpp"
#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "graph.hpp"

namespace CMR {
namespace Sep {

/// Exact primal blossom separation as per Letchford and Lodi's algorithm.
class ExBlossoms {
public:
    ExBlossoms(const std::vector<Graph::Edge> &_edges,
               Data::BestGroup &b_dat, Data::SupportGroup &s_dat,
               CutQueue<ex_blossom> &_blossom_q);
    
    bool find_cuts();
    
private:
    const std::vector<Graph::Edge> &edges;
    const Data::BestGroup &best_data;
    Data::SupportGroup &supp_data;
    CutQueue<ex_blossom> &blossom_q;
    
};

}
}




#endif
