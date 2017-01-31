#include "pricer.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
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
using f64 = util::Fixed64;

namespace Eps = Epsilon;

namespace Price {

using d_PrEdge = PrEdge<double>;

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
    bool silent = true;
    
    edge_hash.clear();

    try {
        reg_duals =
        util::make_unique<LP::DualGroup<double>>(false, core_lp,
                                                 core_lp.external_cuts());
    } CMR_CATCH_PRINT_THROW("populating clique pi", err);

    vector<double> &node_pi_est = reg_duals->node_pi_est;
    CCtsp_edgegenerator *current_eg;

    if (piv_stat == PivType::FathomedTour){
        current_eg = &eg_full;
        if (!silent)
            cout << "\tRunning full eg\n";
        cout << "\n\n";
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

    int finished = 0;
    int outercount = 0;
    int total_added = 0;

    double penalty = 0.0;
    double tourlen = core_lp.get_objval();

    vector<d_PrEdge> price_elist;
    vector<d_PrEdge> edge_q;

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
        price_edges(price_elist, reg_duals);

        for (const d_PrEdge &e : price_elist) {
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
                std::sort(edge_q.begin(), edge_q.end());
                try {
                    vector<Graph::Edge> add_batch = pool_chunk(edge_q);
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
            std::sort(edge_q.begin(), edge_q.end());

            try {
                vector<Graph::Edge> add_batch = pool_chunk(edge_q);
                
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
                reg_duals.reset();
                price_edges(edge_q, reg_duals);
            } CMR_CATCH_PRINT_THROW("getting new duals and re-pricing", err);

            penalty = 0.0;

            edge_q.erase(std::remove_if(edge_q.begin(), edge_q.end(),
                                        [&penalty](const d_PrEdge &e)
                                        -> bool {
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
                for (const d_PrEdge &e : edge_q)
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

// void Pricer::exact_lb()
// {
//     runtime_error err("Problem in Pricer::exact_lb.");
    
//     vector<f64> x_node_pi;
//     vector<f64> x_node_pi_est;
//     vector<f64> x_cut_pi;

//     std::unordered_map<Sep::Clique, f64> x_clq_pi;

//     try {
//         core_lp.external_cuts().get_duals(true, core_lp, x_node_pi,
//                                           x_node_pi_est, x_cut_pi, x_clq_pi);
//     } CMR_CATCH_PRINT_THROW("getting exact duals", err);

//     f64 lb = 0.0;

//     for (f64 &pi : x_node_pi){
//         lb += pi;
// //        lb.add_mult(pi, 2);
//     }

//     cout << "\tNode pi rhs: " << lb << "\n";

//     for (auto i = 0; i < x_cut_pi.size(); ++i) {
//         const Sep::HyperGraph &H = core_lp.external_cuts().get_cuts()[i];
//         if (H.cut_type() == CutType::Non)
//             throw logic_error("Non hypergraph cut in Pricer::exact_lb.");

//         if (H.get_sense() == 'G')
//             lb.add_mult(x_cut_pi[i], H.get_rhs());
//         else
//             lb.add_mult(x_cut_pi[i], -H.get_rhs());
//     }

//     cout << "\tComputed f64 initial rhs: " << lb << "\n";

//     vector<PrEdge<double>64> gen_edges;
//     f64 penalty = 0.0;
    
//     int loop1 = 0;
//     int loop2 = 1;

//     bool finished = false;

//     while (!finished) {
//         try {
//             finished = f64_gen_edges(x_node_pi_est, gen_edges, loop1, loop2);
//             f64_price_edges(gen_edges, x_node_pi, x_node_pi_est, x_cut_pi,
//                             x_clq_pi);
//         } CMR_CATCH_PRINT_THROW("in generating/pricing f64 edges", err);

//         for (auto &e : gen_edges)
//             if (e.redcost < 0.0)
//                 penalty += e.redcost;
//     }
    
//     lb += penalty;
//     cout << "\tPenalty: " << penalty << "\n";
//     cout << "\tAfter subtracting penalty: " << lb << "\n";
// }

vector<Graph::Edge> Pricer::pool_chunk(vector<d_PrEdge> &edge_q)
{
    vector<Graph::Edge> result;

    if (edge_q.size() <= AddBatch) {
        result.reserve(edge_q.size());
        for (const d_PrEdge &e : edge_q)
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

/**
 * A rewrite of the unexported Concorde function big_generate_edges from 
 * ex_price.c. Scans through the edges of the complete graph, building a list 
 * of edges that may have negative reduced cost.
 */
// bool Pricer::f64_gen_edges(const vector<f64> &node_pi_est,
//                            vector<PrEdge<double>64> &gen_edges,
//                            int &loop1, int &loop2)
// {
//     int ncount = inst.node_count();
//     int i = loop1;
//     int j = loop2;
//     int first = 1;

//     gen_edges.clear();

//     if (i >= ncount)
//         return true;

//     for(; i < ncount; ++i) {        
//         int stop = ncount;
//         if (first == 0)
//             j = i + 1;
//         first = 0;
//         for(; j < stop; ++j) {
//             int end = j;
//             f64 rc = inst.edgelen(i, j) - node_pi_est[i] - node_pi_est[j];
//             if (rc < 0.0) {
//                 gen_edges.emplace_back(i, end, rc);
//                 if (gen_edges.size() == f64Batch) {
//                     loop1 = i;
//                     loop2 = j + 1;
//                     return false;
//                 }
//             }
//         }
//     }

//     loop1 = ncount;
//     loop2 = ncount;
//     return true;
// }

// void Pricer::f64_price_edges(vector<PrEdge<double>64> &target_edges,
//                              vector<f64> &node_pi,
//                              vector<f64> &node_pi_est,
//                              vector<f64> &cut_pi,
//                              std::unordered_map<Sep::Clique, f64> &clique_pi)
// {
//     int ncount = inst.node_count();
//     vector<Graph::Edge> temp_elist;
//     vector<PrEdge<double>> tmp_predges;
    
//     for (PrEdge<double>64 &e : target_edges) {
//         e.redcost = inst.edgelen(e.end[0], e.end[1]) - node_pi[e.end[0]]
//         - node_pi[e.end[1]];
//         temp_elist.emplace_back(e.end[0], e.end[1], 0.0);
//         tmp_predges.emplace_back(e.end[0], e.end[1]);
//     }

//     Graph::AdjList price_adjlist(ncount, temp_elist);
//     vector<Graph::Node> &price_nodelist = price_adjlist.nodelist;
//     const vector<int> &def_tour = ext_cuts.get_cbank().ref_tour();
//     int marker = 0;

//     for (const auto &kv : clique_pi) {
//         const Sep::Clique &clq = kv.first;
//         f64 pival = kv.second;

//         if (pival != 0.0) {
//             f64 add_back = pival + pival;
//             ++marker;

//             for (int j : clq.node_list(def_tour)) {
//                 for (Graph::AdjObj &nbr : price_nodelist[j].neighbors)
//                     if (price_nodelist[nbr.other_end].mark == marker)
//                         target_edges[nbr.edge_index].redcost += add_back;
                    
//                 price_nodelist[j].mark = marker;
//             }
//         }
//     }

//     const vector<Sep::HyperGraph> &cutlist = ext_cuts.get_cuts();
//     vector<int> rmatind;
//     vector<double> rmatval;

//     for (int i = 0; i < cutlist.size(); ++i) {
//         f64 pival = cut_pi[i];
//         const Sep::HyperGraph &H = cutlist[i];
        
//         if (H.cut_type() == CutType::Non)
//             throw logic_error("Called pricing with Non HyperGraph present.");
        
//         if (H.cut_type() != CutType::Domino)
//             continue;

//         if (pival == 0)
//             continue;
        
//         try {
//             H.get_coeffs(tmp_predges, rmatind, rmatval);
//         } catch (const exception &e) {
//             cerr << e.what() << "\n";
//             throw runtime_error("Couldn't get price edge domino coeffs.");
//         }

//         for (int j = 0; j < rmatind.size(); ++j)
//             target_edges[rmatind[j]].redcost.add_mult(pival, -rmatval[j]);
//     }
// }


}
}
