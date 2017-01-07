#include "pricer.hpp"
#include "err_util.hpp"



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
    gen_max(EstBatch + ScaleBatch * inst.node_count()),
    gen_elist(vector<int>(2 * gen_max)), gen_elen(gen_max),
    node_pi(vector<double>(inst.node_count())),
    node_pi_est(vector<double>(inst.node_count())),
    cut_pi(vector<double>(relax.num_rows() - inst.node_count()))
{
    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);
    
    if (CCtsp_init_edgegenerator(&eg_inside, inst.node_count(), inst.ptr(),
                                 (CCtsp_genadj *) NULL, Nearest, 1, &rstate))
        throw runtime_error("");

} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Pricer constructor failed.");
}

Pricer::~Pricer() { CCtsp_free_edgegenerator(&eg_inside); }

//before invoking full_scan or partial_scan, this function populates pi info
//for the degree eqns and cuts
//it attempts to reproduce static int pricing_duals from concorde/TSP/tsp_lp.c
ScanStat Pricer::add_edges(PivType piv_stat)
{
    runtime_error err("Problem in Pricer::add_edges");
    ScanStat result = ScanStat::Partial;

    const vector<Sep::HyperGraph> &cutlist = ext_cuts.get_cuts();
    const Sep::CliqueBank &clq_bank = ext_cuts.get_cbank();
    const vector<int> &def_tour = clq_bank.ref_tour();
    
    try { //get updated pi values from the lp
        relax.get_pi(node_pi, 0, inst.node_count() - 1);
        node_pi_est = node_pi;
        
        relax.get_pi(cut_pi, inst.node_count(), relax.num_rows() - 1);
        
        clique_pi.clear();
        clique_pi.reserve(clq_bank.size());
    } CMR_CATCH_PRINT_THROW("clearing/getting/updating pi values", err);

    //get clique_pi for non-domino cuts
    for (int i = 0; i < cutlist.size(); ++i) {
        const Sep::HyperGraph &H = cutlist[i];
        if (H.cut_type() != Sep::HyperGraph::Type::Standard)
            continue;

        double pival = cut_pi[i];

        for (const Sep::Clique::Ptr &clq_ref : H.get_cliques()) {
            if (clique_pi.count(*clq_ref) == 0)
                clique_pi[*clq_ref] = 0.0;
            
            clique_pi[*clq_ref] += pival;
        }
    }

    //use clique_pi to build node_pi for all cliques with nonzero pi
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

    //now get node_pi_est for domino cuts, skipping standard ones
    for (int i = 0; i < cutlist.size(); ++i) {
        const Sep::HyperGraph &H = cutlist[i];
        if (H.cut_type() != Sep::HyperGraph::Type::Domino)
            continue;

        double pival = cut_pi[i];
        if (pival <= 0.0)
            continue;
        
        const Sep::Clique::Ptr &handle_ref = H.get_cliques()[0];
        for (const Segment &seg : handle_ref->seg_list()) {
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                node_pi_est[node] += pival;
            }
        }

        for (const Sep::Tooth::Ptr &T : H.get_teeth())
            for (const Sep::Clique &tpart : T->set_pair())
                for (const Segment &seg : tpart.seg_list())
                    for (int k = seg.start; k <= seg.end; ++k) {
                        int node = def_tour[k];

                        node_pi_est[node] += pival;
                    }
    }

    //now run appropriate pricing scan

    return result;    
}

ScanStat Pricer::partial_scan()
{
    runtime_error err("Problem in Pricer::partial_scan");
    ScanStat result = ScanStat::Partial;

    if (CCtsp_reset_edgegenerator(&eg_inside, &node_pi_est[0], 1))
        throw err;

    int num_gen = 0;
    int finished = 0;
    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);

    if (CCtsp_generate_edges(&eg_inside, gen_max, &num_gen,
                             &gen_elist[0], &gen_elen[0], &finished, 1,
                             &rstate))
        throw err;

    //need to add access to current graph: scroll thru gen_elist seeing if
    //edges found are not in graph yet.
    //probably a good time to get rid of the edgehash from Graph.
    //then call price_list-esque thing on the edges found.
    

    return result;
}

ScanStat Pricer::full_scan()
{
    runtime_error err("Problem in Pricer::full_scan");
    ScanStat result = ScanStat::Full;

    CCtsp_edgegenerator eg;
    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);

    if (CCtsp_init_edgegenerator(&eg, inst.node_count(), inst.ptr(),
                                 (CCtsp_genadj *) NULL,
                                 CCtsp_PRICE_COMPLETE_GRAPH, 1, &rstate)) {
        cerr << "CCtsp_init_edgegenerator failed.\n";
        throw err;
    }
    
    auto cleanup = util::make_guard([&eg] {
        CCtsp_free_edgegenerator(&eg);
    });

    if (CCtsp_reset_edgegenerator(&eg, &node_pi_est[0], 1))
        throw err;

    
    

    return result;
}

}
}
