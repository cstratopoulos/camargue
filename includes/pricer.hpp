#ifndef CMR_PRICER_H
#define CMR_PRICER_H

#include "lp_interface.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "util.hpp"

namespace CMR {

/** Namespace for manners related to pricing sets of edges. */
namespace Price {

constexpr int nearest_factor = 50;

/** Class for pricing edges not in the CoreLP. */
class Pricer {
public:
    Pricer(const LP::Relaxation &_relax, const Data::Instance &_inst,
           const Sep::ExternalCuts &_ext_cuts);

private:
    const LP::Relaxation &relax;
    const Data::Instance &inst;
    const Sep::ExternalCuts &ext_cuts;
    
    int last_ind;

    util::c_array_ptr nearest_elist;
    int nearest_ecount;

};

}
}

#endif
