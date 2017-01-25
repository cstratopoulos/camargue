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
/**
 * @param[in] _relax the Relaxation for grabbing dual values
 * @param[in] _inst the TSP instance for generating edges
 * @param[in] _ext_cuts the HyperGraph representation of the cuts in 
 * \p _relax.
 */
Pricer::Pricer(LP::CoreLP &core, const Data::Instance &_inst,
               Data::GraphGroup &graphgroup) try :
    core_lp(core), inst(_inst), ext_cuts(core.external_cuts()),
    graph_group(graphgroup),
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
        throw runtime_error("CCtsp_init_edgegenerator(50) failed.");

    if (CCtsp_init_edgegenerator(&eg_full, inst.node_count(), inst.ptr(),
                                 (CCtsp_genadj *) NULL,
                                 CCtsp_PRICE_COMPLETE_GRAPH, 1, &rstate))
        throw runtime_error("CCtsp_init_edgenerator(complete) failed.");

} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Pricer constructor failed.");
}

Pricer::~Pricer()
{
    CCtsp_free_edgegenerator(&eg_inside);
    CCtsp_free_edgegenerator(&eg_full);
}


/**
 * @param[in] piv_stat the PivType from the solution when edge generation 
 * is called. The partial edge set is scanned iff \p piv_stat is 
 * PivType::Tour. Full edge set may be scanned if \p piv_stat is 
 * PivType::FathomedTour, or PivType::Tour with a small enough graph, or
 * if a partial scan finds no edges to add. 
 * @returns a ScanStat indicating the outcome of the pricing scan.
 */
ScanStat Pricer::gen_edges(LP::PivType piv_stat)
{
    runtime_error err("Problem in Pricer::gen_edges");
    ScanStat result =
    (piv_stat == LP::PivType::Tour) ? ScanStat::PartOpt : ScanStat::FullOpt;
    bool silent = false;
    
    edge_hash.clear();

    try {
        ext_cuts.get_duals(core_lp, node_pi, node_pi_est, cut_pi, clique_pi);
    } CMR_CATCH_PRINT_THROW("populating clique pi", err);

    CCtsp_edgegenerator *current_eg;

    if (piv_stat == PivType::FathomedTour){
        current_eg = &eg_full;
        if (!silent)
            cout << "\tRunning full eg\n";
    } else if (piv_stat == PivType::Tour){
        current_eg = &eg_inside;
        if (!silent)
            cout << "\tRunning inside eg\n";
    } else
        throw logic_error("Tried to run pricing on non tour.");

    if (CCtsp_reset_edgegenerator(current_eg, &node_pi_est[0], silent)) {
        cerr << "CCtsp_reset_edgegenerator failed.\n";
        throw err;
    }

    Graph::CoreGraph &core_graph = graph_group.core_graph;
    CCrandstate rstate;

    CCutil_sprand(inst.seed(), &rstate);

    edge_q.clear();

    int finished = 0;
    int outercount = 0;
    int total_added = 0;

    double penalty = 0.0;
    double tourlen = core_lp.get_objval();

    while (!finished) {
        if (!silent)
            cout << "\tEntering EG loop, pass " << ++outercount << "\n\t";

        int num_gen = 0;

        price_elist.clear();

        if (CCtsp_generate_edges(current_eg, gen_max, &num_gen, &gen_elist[0],
                                 &gen_elen[0], &finished, silent, &rstate)) {
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

        if (!silent)
            cout << "\t Pricing candidates\n";
        price_edges(price_elist, false);

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

        if (!silent)
            cout << "\t" << edge_q.size() << " edges in edge_q\n"
                 << "\t" << penalty << " penalty accrued, finished: "
                 << finished << "\n";

        if (piv_stat == PivType::Tour) {
            if (!edge_q.empty()) {
                sort_q();
                try {
                    vector<Graph::Edge> add_batch = get_pool_chunk();
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
        double new_objval = 0.0;
        int num_added = 0;
        int innercount = 0;
        
        while (((!finished && edge_q.size() >= PoolSize) ||
                (finished && penalty < -MaxPenalty && !edge_q.empty()))) {
            if (!silent)
                cout << "\t\tOpt tour inner price loop, pass " << ++innercount
                     << "\n";
            sort_q();

            try {
                vector<Graph::Edge> add_batch = get_pool_chunk();
                
                num_added = add_batch.size();
                total_added += num_added;
                core_lp.add_edges(add_batch);
                core_lp.primal_opt();
                new_objval = core_lp.get_objval();
            } CMR_CATCH_PRINT_THROW("adding edges to lp and optimizing", err);

            if (!silent)
                cout << "\t\tAdded " << num_added
                     << " edges in opt solution.\n";

            if (std::abs(new_objval - tourlen) >= Eps::Zero) {
                if (innercount == 1)
                    if (!silent)
                        cout << "\tTour no longer optimal after adding edges\n";
                result = ScanStat::Full;
            } else {
                if (!silent)
                    cout << "\t\tTour still optimal after adding edges.\n";
                result = ScanStat::FullOpt;
            }

            try {
                price_edges(edge_q, true);
            } CMR_CATCH_PRINT_THROW("getting new duals and re-pricing", err);

            penalty = 0.0;

            edge_q.erase(std::remove_if(edge_q.begin(), edge_q.end(),
                                        [&penalty](const PrEdge &e) -> bool {
                                            if (e.redcost < 0.0)
                                                penalty += e.redcost;
                                            return e.redcost > - Eps::Zero;
                                        }),
                         edge_q.end());


            if (num_added > 0) {
                if (CCtsp_reset_edgegenerator(current_eg, &node_pi_est[0],
                                              silent)) {
                    cerr << "CCtsp_reset_edgegenerator failed.\n";
                    throw err;
                }
                
                finished = 0;
                edge_hash.clear();
                for (const PrEdge &e : edge_q)
                    edge_hash.add(e.end[0], e.end[1], 1);
            }

            if (!silent)
                cout << "\t\tInner penalty: "
                     << penalty << ", new edge_q size "
                     << edge_q.size() << ", finished: " << finished << "\n";
        }
    }
    
    cout << "\tAdded " << total_added << " edges, " << result << "\n";
    return result;
}

vector<Graph::Edge> Pricer::get_pool_chunk()
{
    vector<Graph::Edge> result;

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

void Pricer::price_edges(vector<PrEdge> &target_edges, bool compute_duals)
{
    runtime_error err("Problem in Pricer::price_edges.");
    
    if (compute_duals)
        try {
            ext_cuts.get_duals(core_lp, node_pi, node_pi_est, cut_pi,
                               clique_pi);
        } CMR_CATCH_PRINT_THROW("Couldn't get duals.", err);

    for (PrEdge &e : target_edges)
        e.redcost = inst.edgelen(e.end[0], e.end[1]) - node_pi[e.end[0]]
        - node_pi[e.end[1]];
    
    Graph::AdjList price_adjlist;

    try  {
        price_adjlist = Graph::AdjList(inst.node_count(), target_edges);
    } CMR_CATCH_PRINT_THROW("Couldn't build price adjlist.", err);

    vector<Graph::Node> &price_nodelist = price_adjlist.nodelist;
    
    const std::vector<int> &def_tour = ext_cuts.get_cbank().ref_tour();
    int marker = 0;

    
    for (const std::pair<Sep::Clique, double> &kv : clique_pi) {
        const Sep::Clique &clq = kv.first;
        double pival = kv.second;

        if (pival != 0.0) {
            double add_back = 2 * pival;
            ++marker;
            
            for (int j : clq.node_list(def_tour)) {
                for (Graph::AdjObj &nbr : price_nodelist[j].neighbors)
                    if (price_nodelist[nbr.other_end].mark == marker)
                        target_edges[nbr.edge_index].redcost += add_back;
                    
                price_nodelist[j].mark = marker;
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
        if (H.cut_type() != CutType::Domino)
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
