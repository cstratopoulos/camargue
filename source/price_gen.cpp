/**
 * @file
 * @brief Pricer implementations for generating negative reduced cost edges.
 */

#include "pricer.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <utility>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

using std::vector;

using std::cout;
using std::cerr;
using std::endl;

using std::exception;
using std::runtime_error;

namespace CMR {

using LP::PivType;
using CutType = Sep::HyperGraph::Type;

namespace Eps = Epsilon;

namespace Price {

using d_PrEdge = PrEdge<double>;

struct Pricer::edgegen_impl {
    edgegen_impl(const Data::Instance &inst, int price_mode);
    ~edgegen_impl();

    void gen_reset(vector<double> &node_pi_est);
    void gen_edges(int max_ngen, int &num_gen, vector<int> &generated_elist,
                   vector<int> &generated_elen, int &finished,
                   CCrandstate &rstate);

    CCtsp_edgegenerator eg;
};

Pricer::edgegen_impl::edgegen_impl(const Data::Instance &inst, int price_mode)
{
    CCrandstate rstate;
    CCutil_sprand(inst.seed(), &rstate);

    if (CCtsp_init_edgegenerator(&eg, inst.node_count(), inst.ptr(), NULL,
                                 price_mode, 1, &rstate))
        throw runtime_error("CCtsp_init_edgegenerator failed");
}

Pricer::edgegen_impl::~edgegen_impl() { CCtsp_free_edgegenerator(&eg); }

void Pricer::edgegen_impl::gen_reset(vector<double> &node_pi_est)
{
    if (CCtsp_reset_edgegenerator(&eg, &node_pi_est[0], 1))
        throw runtime_error("CCtsp_reset_edgegenerator failed");
}

void Pricer::edgegen_impl::gen_edges(int max_ngen, int &num_gen,
                                     vector<int> &generated_elist,
                                     vector<int> &generated_elen,
                                     int &finished, CCrandstate &rstate)
{
    if (CCtsp_generate_edges(&eg, max_ngen, &num_gen,
                             &generated_elist[0], &generated_elen[0],
                             &finished, 1, &rstate))
        throw runtime_error("CCtsp_generate_edges failed");
}

/**
 * @param[in] _relax the Relaxation for grabbing dual values
 * @param[in] _inst the TSP instance for generating edges
 * @param[in] _ext_cuts the HyperGraph representation of the cuts in
 * \p _relax.
 */
Pricer::Pricer(LP::CoreLP &core, const Data::Instance &_inst,
               Graph::CoreGraph &core_graph_) try :
    core_lp(core), inst(_inst), ext_cuts(core.external_cuts()),
    core_graph(core_graph_),
    gen_max(EstBatch + ScaleBatch * inst.node_count()),
    gen_elist(vector<int>(2 * gen_max)), gen_elen(gen_max),
    eg_inside(util::make_unique<edgegen_impl>(_inst, 50)),
    eg_full(util::make_unique<edgegen_impl>(_inst, CCtsp_PRICE_COMPLETE_GRAPH)),
    edge_hash(gen_max)
{} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Pricer constructor failed.");
}

Pricer::~Pricer() { }


/**
 * @param[in] piv_stat the PivType from the solution when edge generation
 * is called. The partial edge set is scanned iff \p piv_stat is
 * PivType::Tour. Full edge set may be scanned if \p piv_stat is
 * PivType::FathomedTour, or PivType::Tour with a small enough graph, or
 * if a partial scan finds no edges to add.
 * @returns a ScanStat indicating the outcome of the pricing scan.
 */
ScanStat Pricer::gen_edges(LP::PivType piv_stat, bool try_elim)
{
    runtime_error err("Problem in Pricer::gen_edges");
    ScanStat result =
    (piv_stat == LP::PivType::Tour) ? ScanStat::PartOpt : ScanStat::FullOpt;

    edge_hash.clear();

    try { util::ptr_reset(reg_duals, false, core_lp, ext_cuts); }
    CMR_CATCH_PRINT_THROW("populating clique pi", err);

    edgegen_impl *current_eg;

    if (piv_stat == PivType::FathomedTour){
        current_eg = eg_full.get();
        if (verbose > 1)
            cout << "\tRunning full eg\n";
    } else if (piv_stat == PivType::Tour){
        current_eg = eg_inside.get();
        if (verbose > 1)
            cout << "\tRunning inside eg\n";
    } else
        throw runtime_error("Tried to run pricing on non tour.");

    try { current_eg->gen_reset(reg_duals->node_pi_est); }
    CMR_CATCH_PRINT_THROW("resetting generator", err);

    CCrandstate rstate;

    CCutil_sprand(inst.seed(), &rstate);

    int finished = 0;
    int outercount = 0;
    int inner_total = 0;
    int total_added = 0;

    double penalty = 0.0;

    double upper_bound = core_lp.global_ub();
    //double lower_bound = core_lp.get_objval();

    vector<d_PrEdge> price_elist;
    vector<d_PrEdge> edge_q;

    while (!finished) {
        ++outercount;

        int num_gen = 0;

        price_elist.clear();

        try {
            current_eg->gen_edges(gen_max, num_gen, gen_elist, gen_elen,
                                  finished, rstate);
        } CMR_CATCH_PRINT_THROW("generating edges", err);

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

        price_edges(price_elist, reg_duals, true);

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

        if (piv_stat == PivType::Tour) {
            if (!edge_q.empty()) {
                std::sort(edge_q.begin(), edge_q.end());
                try {
                    vector<Graph::Edge> add_batch = pool_chunk(edge_q);
                    core_lp.add_edges(add_batch, true);
                } CMR_CATCH_PRINT_THROW("adding edges for aug tour", err);

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

        while ((!finished && edge_q.size() >= PoolSize) ||
               (finished && penalty < -MaxPenalty && !edge_q.empty())) {
            ++inner_total;
            std::sort(edge_q.begin(), edge_q.end());

            try {
                vector<Graph::Edge> add_batch = pool_chunk(edge_q);

                num_added = add_batch.size();
                total_added += num_added;
                core_lp.add_edges(add_batch, true);
                core_lp.primal_opt();
                new_objval = core_lp.get_objval();
            } CMR_CATCH_PRINT_THROW("adding edges to lp and optimizing", err);

            if (new_objval >= upper_bound - 0.9)
                result = ScanStat::FullOpt;
            else
                result = ScanStat::Full;

            try {
                reg_duals.reset();
                price_edges(edge_q, reg_duals, true);
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
                try { current_eg->gen_reset(reg_duals->node_pi_est); }
                CMR_CATCH_PRINT_THROW("resetting generator", err);

                finished = 0;
                edge_hash.clear();
                try {
                    for (const d_PrEdge &e : edge_q)
                        edge_hash.add(e.end[0], e.end[1], 1);
                } CMR_CATCH_PRINT_THROW("adding q edges to hash", err);
            }
        }
    }

    if (verbose > 1)
        cout << "Added " << total_added << " edges in "
             << outercount << " passes, " << inner_total << " inner loops"
             << endl;

    if (verbose)
        cout << "LP opt objval: " << core_lp.get_objval() << ", "
             << core_lp.num_cols() << " columns" << endl;

    if (try_elim && result != ScanStat::FullOpt && core_lp.dual_feas()) {
        try { elim_edges(false); }
        CMR_CATCH_PRINT_THROW("eliminating after gen", err);
    }

    return result;
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


}
}
