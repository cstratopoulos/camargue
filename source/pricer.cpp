#include "pricer.hpp"
#include "err_util.hpp"

extern "C" {
#include <concorde/INCLUDE/edgegen.h>
}

#include <iostream>
#include <stdexcept>
#include <utility>

using std::vector;

using std::cout;
using std::cerr;

using std::exception;
using std::runtime_error;
using std::logic_error;

namespace CMR {

using LP::PivType;

namespace Price {

Pricer::Pricer(const LP::Relaxation &_relax, const Data::Instance &_inst,
               const Sep::ExternalCuts &_ext_cuts) try :
    relax(_relax), inst(_inst), ext_cuts(_ext_cuts),
    node_pi(vector<double>(inst.node_count())),
    cut_pi(vector<double>(relax.num_rows() - inst.node_count())),
    last_ind(0)
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

ScanStat Pricer::add_edges(PivType piv_stat)
{
    runtime_error err("Problem in Pricer::add_edges");
    ScanStat result = ScanStat::Partial;

    const vector<Sep::HyperGraph> &cutlist = ext_cuts.get_cuts();
    const Sep::CliqueBank &clq_bank = ext_cuts.get_cbank();
    const vector<int> &def_tour = clq_bank.ref_tour();
    
    try {
        relax.get_pi(node_pi, 0, inst.node_count() - 1);
        node_pi_est = node_pi;
        
        relax.get_pi(cut_pi, inst.node_count(), relax.num_rows() - 1);
        
        clique_pi.clear();
        clique_pi.reserve(clq_bank.size());
    } CMR_CATCH_PRINT_THROW("clearing/getting/updating pi values", err);

    for (int i = 0; i < cutlist.size(); ++i) {
        const Sep::HyperGraph &H = cutlist[i];
        if (H.cut_type() == Sep::HyperGraph::Type::Domino)
            continue;

        double pival = cut_pi[i];

        for (const Sep::Clique::Ptr &clq_ref : H.get_cliques()) {
            if (clique_pi.count(*clq_ref) == 0)
                clique_pi[*clq_ref] = 0.0;
            
            clique_pi[*clq_ref] += pival;
        }
    }

    for (const std::pair<Sep::Clique, double> &kv : clique_pi) {
        const Sep::Clique &clq = kv.first;
        double pival = kv.second;

        if (pival > 0.0) {
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];
                    
                    node_pi[node] += pival;
                    node_pi_est[node] += pival;
                }
            }            
        } else if (pival < 0.0) {
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];
                    
                    node_pi[node] += pival;
                }
            }
        }
    }

    

    return result;    
}

}
}
