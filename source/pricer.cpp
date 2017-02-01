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
        reg_duals = util::make_unique<LP::DualGroup<double>>(false, core_lp,
                                                             ext_cuts);
    } CMR_CATCH_PRINT_THROW("populating clique pi", err);

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

    if (CCtsp_reset_edgegenerator(current_eg,
                                  &(reg_duals->node_pi_est[0]), silent)) {
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
                if (CCtsp_reset_edgegenerator(current_eg,
                                              &(reg_duals->node_pi_est[0]),
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

f64 Pricer::exact_lb()
{
    runtime_error err("Problem in Pricer::exact_lb");

    cout << "\tObj val: " << core_lp.get_objval() << ", dual feas: "
         << core_lp.dual_feas() << ", getting duals....\n";

    try {
        ex_duals = util::make_unique<LP::DualGroup<f64>>(true, core_lp,
                                                         ext_cuts);
    } CMR_CATCH_PRINT_THROW("constructing exact DualGroup", err);

    vector<f64> &node_pi = ex_duals->node_pi;
    vector<f64> &cut_pi = ex_duals->cut_pi;


    cout << "\tnode_pi size: " << node_pi.size() << "\n";
    cout << "\tcut_pi size: " << cut_pi.size() << "\n";

    const vector<Sep::HyperGraph> &cuts = ext_cuts.get_cuts();

    f64 node_sum = 0.0;
    f64 cut_sum = 0.0;

    for (const f64 &pi : ex_duals->node_pi) {
        util::add_mult(node_sum, pi, 2);
    }
    cout << "\tnode pi sum:\t" << node_sum << "\n";
    
    for (int i = 0; i < cuts.size(); ++i) {
        const Sep::HyperGraph &H = cuts[i];
        if (H.cut_type() == CutType::Non)
            throw logic_error("Tried to get exact_lb with non cut present.");

        cout << "\tAdding mult " << cut_pi[i] << ", " << H.get_rhs() << "\n";
        util::add_mult(cut_sum, cut_pi[i], H.get_rhs());
    }

    cout << "\tcut pi sum:\t" << cut_sum << "\n";

    vector<PrEdge<f64>> graph_edges;

    for (const Graph::Edge &e : graph_group.core_graph.get_edges())
        graph_edges.emplace_back(e.end[0], e.end[1]);

    price_edges(graph_edges, ex_duals);

    f64 rc_sum = 0.0;
    
    for (const PrEdge<f64> &e : graph_edges)
        if (e.redcost < 0.0)
            rc_sum -= e.redcost;

    cout << "\tred cost sum:\t" << rc_sum << "\n";

    f64 bound = node_sum + cut_sum - rc_sum;

    cout << "\tFinal bound:\t" << bound << "\n";    
    return bound;

    // vector<PrEdge<f64>> gen_edges;
    // f64 penalty = 0.0;

    // bool finished = false;
    // int loop1 = 0;
    // int loop2 = 1;

    // while (!finished) {
    //     try {
    //         finished = scan_edges(gen_edges, loop1, loop2);
    //         price_edges(gen_edges, ex_duals);
    //     } CMR_CATCH_PRINT_THROW("generating/pricing edges", err);

    //     for (const PrEdge<f64> &e : gen_edges)
    //         if (e.redcost < 0.0)
    //             penalty += e.redcost;
    // }

    // cout << "\tAccrued penalty of " << penalty << "\n";
    // bound = rhs_sum + penalty;
    // cout << "\tBound: " << bound << "\n";
    // return bound;
}


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

bool Pricer::scan_edges(vector<PrEdge<f64>> &gen_edges, int &loop1,
                        int &loop2)
{
    if (!ex_duals)
        throw logic_error("Tried to scan edges without exact duals.");

    const vector<f64> &node_pi_est = ex_duals->node_pi_est;
    const Graph::AdjList &alist = graph_group.core_graph.get_adj();

    int ncount = inst.node_count();
    int i = loop1;
    int j = loop2;
    int first = 1;

    gen_edges.clear();

    if (i >= ncount)
        return true;

    for (; i < ncount; ++i) {
        for (const Graph::AdjObj &a : alist.nodelist[i].neighbors) {
            int j = a.other_end;
            if (j > i) {
                f64 len = a.val;
                f64 rc = len - node_pi_est[i] - node_pi_est[j];
                if (rc < 0.0)
                    gen_edges.emplace_back(i, j, rc);
                if (gen_edges.size() == f64Batch) {
                    loop1 = i;
                    return false;
                }
            }
        }
    }

    // for(; i < ncount; ++i) {        
    //     int stop = ncount;
    //     if (first == 0)
    //         j = i + 1;
    //     first = 0;
    //     for(; j < stop; ++j) {
    //         int end = j;
    //         f64 rc = inst.edgelen(i, j) - node_pi_est[i] - node_pi_est[j];
    //         if (rc < 0.0) {
    //             gen_edges.emplace_back(i, end, rc);
    //             if (gen_edges.size() == f64Batch) {
    //                 loop1 = i;
    //                 loop2 = j + 1;
    //                 return false;
    //             }
    //         }
    //     }
    // }

    loop1 = ncount;
    loop2 = ncount;
    return true;    
}


}
}
