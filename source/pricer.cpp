#include "pricer.hpp"
#include "err_util.hpp"

#include <algorithm>
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
using CutType = Sep::HyperGraph::Type;

namespace Eps = Epsilon;

namespace Price {

Pricer::Pricer(LP::CoreLP &core, const Data::Instance &_inst,
               const Sep::ExternalCuts &_ext_cuts,
               Data::GraphGroup &graphgroup) try :
    core_lp(core), inst(_inst), ext_cuts(_ext_cuts), graph_group(graphgroup),
    gen_max(EstBatch + ScaleBatch * inst.node_count()),
    gen_elist(vector<int>(2 * gen_max)), gen_elen(gen_max),
    node_pi(vector<double>(inst.node_count())),
    node_pi_est(vector<double>(inst.node_count())),
    cut_pi(vector<double>(core_lp.num_rows() - inst.node_count())),
    edge_hash(gen_max)
{
    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);
    
    if (CCtsp_init_edgegenerator(&eg_inside, inst.node_count(), inst.ptr(),
                                 (CCtsp_genadj *) NULL, Nearest, 1, &rstate))
        throw runtime_error("");

    if (CCtsp_init_edgegenerator(&eg_full, inst.node_count(), inst.ptr(),
                                 (CCtsp_genadj *) NULL,
                                 CCtsp_PRICE_COMPLETE_GRAPH, 1, &rstate))
        throw runtime_error("");

} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Pricer constructor failed.");
}

Pricer::~Pricer()
{
    CCtsp_free_edgegenerator(&eg_inside);
    CCtsp_free_edgegenerator(&eg_full);
}

ScanStat Pricer::gen_edges(LP::PivType piv_stat)
{
    runtime_error err("Problem in Pricer::gen_edges");
    
    edge_hash.clear();

    try {
        get_duals();
    } CMR_CATCH_PRINT_THROW("populating clique pi", err);

    CCtsp_edgegenerator *current_eg;

    if (piv_stat == PivType::FathomedTour){
        current_eg = &eg_full;
        cout << "\tRunning full eg\n";
    } else if (piv_stat == PivType::Tour){
        current_eg = &eg_inside;
        cout << "\tRunning inside eg\n";
    } else
        throw logic_error("Tried to run pricing on non tour.");

    if (CCtsp_reset_edgegenerator(current_eg, &node_pi_est[0], 0)) {
        cerr << "CCtsp_reset_edgegenerator failed.\n";
        throw err;
    }

    GraphUtils::CoreGraph &core_graph = graph_group.core_graph;
    CCrandstate rstate;

    CCutil_sprand(inst.seed(), &rstate);

    edge_q.clear();

    int finished = 0;
    int outercount = 0;

    double penalty = 0.0;

    while (!finished) {
        cout << "\tEntering EG loop, pass " << ++outercount << "\n\t";

        int num_gen = 0;

        price_elist.clear();

        if (CCtsp_generate_edges(current_eg, gen_max, &num_gen, &gen_elist[0],
                                 &gen_elen[0], &finished, 0, &rstate)) {
            cerr << "CCtsp_generate_edges failed.\n";
            throw err;
        }

        for (int i = 0; i < num_gen; ++i) {
            int e0 = gen_elist[2 * i];
            int e1 = gen_elist[(2 * i) + 1];
            int hashval = edge_hash.get_val(e0, e1);
            if (hashval == -1 && core_graph.find_edge_ind(e0, e1) == -1) {
                try {
                    price_elist.emplace_back(e0, e1);
                    edge_hash.add(e0, e1, 1);
                } CMR_CATCH_PRINT_THROW("adding new canididate edge", err);
            }
        }

        cout << "\t Pricing candidates\n";
        price_candidates();

        for (const PrEdge &e : price_elist) {
            if (e.redcost < 0.0)
                penalty += e.redcost;
            if (e.redcost <= -Eps::Zero) {
                try {
                    edge_q.push_back(e);
                } CMR_CATCH_PRINT_THROW("enqueueing edge for addition", err);
            } else {
                edge_hash.erase(e.end[0], e.end[1]);
            }
        }

        cout << "\t" << edge_q.size() << " edges in edge_q\n"
             << "\t" << penalty << " penalty accrued, finished: "
             << finished << "\n";

        if (piv_stat == PivType::Tour) {
            if (!edge_q.empty()) {
                sort_q();
                try {
                    vector<Edge> add_batch = get_pool_chunk();
                    core_lp.add_edges(add_batch);
                } CMR_CATCH_PRINT_THROW("adding edges for aug tour", err);
                cout << "\tFound and added edges for aug tour.\n";
                
                return ScanStat::Partial;                
            } else if (finished) {
                return ScanStat::PartOpt;
            } else
                continue;
        }

        //this loop is only entered for fathomed tours
        //we use the concorde termination criteria, but also terminate within
        //the loop if adding edges and optimizing takes us to a new lp solution
        double opt_tourlen = core_lp.get_objval();
        double new_objval = opt_tourlen;
        int num_added = 0;
        int innercount = 0;
        
        while (((!finished && edge_q.size() >= PoolSize) ||
                (finished && penalty < -MaxPenalty && !edge_q.empty()))) {
            cout << "\t\tOpt tour inner price loop, pass " << ++innercount
                 << "\n";
            sort_q();

            try {
                vector<Edge> add_batch = get_pool_chunk();
                
                num_added = add_batch.size();
                core_lp.add_edges(add_batch);
                new_objval = core_lp.opt_objval();
            } CMR_CATCH_PRINT_THROW("adding edges to lp and optimizing", err);

            cout << "\t\tAdded " << num_added << " edges in opt solution.\n";

            if (std::abs(new_objval - opt_tourlen) >= Eps::Zero) {
                cout << "\tTour no longer optimal after adding edges.\n";
                return ScanStat::Full;
            }

            cout << "\t\tTour still optimal after adding edges.\n";

            try {
                get_duals();
                price_candidates();
            } CMR_CATCH_PRINT_THROW("getting new duals and re-pricing", err);

            penalty = 0.0;

            edge_q.erase(std::remove_if(edge_q.begin(), edge_q.end(),
                                        [&penalty](const PrEdge &e) -> bool {
                                            if (e.redcost < 0.0)
                                                penalty += e.redcost;
                                            return e.redcost > - Eps::Zero;
                                        }),
                         edge_q.end());

            cout << "\t\tInner penalty: " << penalty << ", new edge_q size "
                 << edge_q.size() << "\n";

            if (num_added > 0) {
                if (CCtsp_reset_edgegenerator(current_eg, &node_pi_est[0],
                                              0)) {
                    cerr << "CCtsp_reset_edgegenerator failed.\n";
                    throw err;
                }
                
                finished = 0;
                edge_hash.clear();
                for (const PrEdge &e : edge_q)
                    edge_hash.add(e.end[0], e.end[1], 1);
            }        
        }
    }
    cout << "\t Finished: " << finished << ", tour seems genuinely optimal.\n";
    return ScanStat::FullOpt;
}

vector<Edge> Pricer::get_pool_chunk()
{
    vector<Edge> result;

    if (edge_q.size() <= AddBatch) {
        result.reserve(edge_q.size());
        for (const PrEdge &e : edge_q)
            result.emplace_back(e.end[0], e.end[1],
                                inst.edgelen(e.end[0], e.end[1]));
        edge_q.clear();
        return result;
    }

    result.reserve(AddBatch);
    for (auto it = edge_q.end() - AddBatch; it != edge_q.end(); ++it)
        result.emplace_back(it->end[0], it->end[1],
                            inst.edgelen(it->end[0], it->end[1]));
    edge_q.resize(edge_q.size() - AddBatch);
    return result;
}

void Pricer::get_duals()
{
    runtime_error err("Problem in Pricer::get_duals");
    
    const vector<Sep::HyperGraph> &cutlist = ext_cuts.get_cuts();
    const Sep::CliqueBank &clq_bank = ext_cuts.get_cbank();
    const vector<int> &def_tour = clq_bank.ref_tour();

    clique_pi.clear();
    
    try { //get updated pi values from the lp
        clique_pi.reserve(clq_bank.size());

        core_lp.get_pi(node_pi, 0, inst.node_count() - 1);
        node_pi_est = node_pi;

        if (!cutlist.empty())
            core_lp.get_pi(cut_pi, inst.node_count(), core_lp.num_rows() - 1);
    } CMR_CATCH_PRINT_THROW("clearing/getting/updating pi values", err);

    //get clique_pi for non-domino cuts
    for (int i = 0; i < cutlist.size(); ++i) {
        const Sep::HyperGraph &H = cutlist[i];
        if (H.cut_type() != CutType::Standard)
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
        if (H.cut_type() != CutType::Domino)
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
}

void Pricer::price_candidates()
{
    if (price_elist.empty() && edge_q.empty())
        return;

    vector<PrEdge> &target_edges = edge_q.empty() ? price_elist : edge_q;
    
    for (PrEdge &e : target_edges)
        e.redcost = inst.edgelen(e.end[0], e.end[1]) - node_pi[e.end[0]]
        - node_pi[e.end[1]];

    GraphUtils::AdjList price_adjlist(inst.node_count(), target_edges);
    vector<GraphUtils::Node> &price_nodelist = price_adjlist.nodelist;
    
    const std::vector<int> &def_tour = ext_cuts.get_cbank().ref_tour();
    int marker = 0;
    
    for (const std::pair<Sep::Clique, double> &kv : clique_pi) {
        const Sep::Clique &clq = kv.first;
        double pival = kv.second;

        if (pival != 0.0) {
            double add_back = 2 * pival;
            ++marker;
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];

                    for (GraphUtils::AdjObj &a :
                         price_nodelist[node].neighbors) {
                        if (price_nodelist[a.other_end].mark == marker)
                            target_edges[a.edge_index].redcost += add_back;
                    }
                    price_nodelist[node].mark = marker;
                }
            }
        }
    }

    const vector<Sep::HyperGraph> &cutlist = ext_cuts.get_cuts();
    vector<int> rmatind;
    vector<double> rmatval;

    for (int i = 0; i < cutlist.size(); ++i) {
        double pival = cut_pi[i];
        if (pival <= 0.0)
            continue;

        const Sep::HyperGraph &H = cutlist[i];
        if (H.cut_type() == CutType::Standard)
            continue;

        try {
            H.get_coeffs(target_edges, rmatind, rmatval);
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            throw runtime_error("Couldn't get price edge domino coeffs.");
        }

        for (int j = 0; j < rmatind.size(); ++j)
            target_edges[rmatind[j]].redcost -= pival * rmatval[j];
    }
}

void Pricer::sort_q()
{
    std::sort(edge_q.begin(), edge_q.end(),
              [](const PrEdge &e, const PrEdge &f) {return f < e; });
}


}
}
