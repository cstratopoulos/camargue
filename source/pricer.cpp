#include "pricer.hpp"

extern "C" {
#include <concorde/INCLUDE/edgegen.h>
}

#include <stdexcept>

using std::exception;
using std::runtime_error;
using std::logic_error;

namespace CMR {
namespace Price {

Pricer::Pricer(const LP::Relaxation &_relax, const Data::Instance &_inst,
               const Sep::ExternalCuts &_ext_cuts) try :
    relax(_relax), inst(_inst), ext_cuts(_ext_cuts), last_ind(0)
{
    int ncount = inst.node_count();
    CCedgegengroup plan;
    CCrandstate rstate;

    CCutil_sprand(inst.seed(), &rstate);
    CCedgegen_init_edgegengroup(&plan);

    plan.nearest = nearest_factor;

    int *elist = (int *) NULL;
    int ecount = 0;

    if (CCedgegen_edges(&plan, ncount, inst.ptr(), NULL, &ecount, &elist,
                        1, &rstate))
        throw runtime_error("Problem in CCedgegen_edges");

    nearest_elist.reset(elist);

    if (ecount >= (ncount * (ncount - 1)) / 2) { //if 50-nearest has full edges
        nearest_elist.reset(nullptr);
        nearest_ecount = 0;
    } else
        nearest_ecount = ecount;
    
} catch (const exception &e) {
    throw runtime_error("Pricer constructor failed.");
}

}
}
