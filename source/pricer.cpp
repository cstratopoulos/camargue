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
using std::logic_error;

namespace CMR {

using LP::PivType;
using CutType = Sep::HyperGraph::Type;
using f64 = util::Fixed64;

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

    try {
        reg_duals = util::make_unique<LP::DualGroup<double>>(false, core_lp,
                                                             ext_cuts);
    } CMR_CATCH_PRINT_THROW("populating clique pi", err);

    edgegen_impl *current_eg;

    if (piv_stat == PivType::FathomedTour){
        current_eg = eg_full.get();
        if (verbose)
            cout << "\tRunning full eg\n";
    } else if (piv_stat == PivType::Tour){
        current_eg = eg_inside.get();
        if (verbose)
            cout << "\tRunning inside eg\n";
    } else
        throw logic_error("Tried to run pricing on non tour.");

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

    cout << "Added " << total_added << " edges in "
         << outercount << " passes, " << inner_total << " inner loops, "
         << result << "\n\t"
         << "Opt objval: " << core_lp.get_objval() << ", "
         << core_lp.num_cols() << " columns" << endl;

    if (try_elim && result != ScanStat::FullOpt && core_lp.dual_feas()) {
        try { elim_edges(false); }
        CMR_CATCH_PRINT_THROW("eliminating after gen", err);
    }

    return result;
}

/**
 * @returns true if adding edges was able to recover feasibility, false
 * if the LP is  genuinely infeasible.
 */
bool Pricer::feas_recover()
{
    using SolStat = LP::SolStat;
    SolStat cur_stat = SolStat::Infeas;

    runtime_error err("Problem in Pricer::feas_recover");

    try {
        core_lp.primal_opt();
        cur_stat = core_lp.get_stat();
    } CMR_CATCH_PRINT_THROW("optimizing the lp", err);

    if (cur_stat != SolStat::Infeas) {
        if (verbose)
            cout << "Pricer::feas_recover was called with stat "
                 << cur_stat << ", returning true." << endl;
        return true;
    } else if (verbose)
        cout << "Primal optimized with infeasible LP" << endl;

    try {
        reg_duals = util::make_unique<LP::DualGroup<double>>(false, core_lp,
                                                             ext_cuts);
    } CMR_CATCH_PRINT_THROW("populating clique pi", err);

    edge_hash.clear();

    int start = 0;
    int end1 = 0;
    int end2 = 1;

    bool finished = false;
    int total_added = 0;
    double penalty = 0.0;

    vector<d_PrEdge> price_elist;
    vector<d_PrEdge> edge_q;

    while (!finished) {
        price_elist.clear();

        try {
            finished = feas_gen_edges(price_elist, start, end1, end2);
            if (verbose)
                cout << "feas_gen_edges start " << start << ", end1 "
                     << end1 << ", end2 " << end2 << ", "
                     << price_elist.size() << " in list, finished "
                     << finished << endl;
        } CMR_CATCH_PRINT_THROW("generating edges", err);

        price_edges(price_elist, reg_duals, false);

        for (const d_PrEdge &e : price_elist) {
            if (e.redcost < 0.0)
                penalty += e.redcost;
            if (e.redcost < -RecoverMaxPen) {
                try {
                    edge_q.push_back(e);
                } CMR_CATCH_PRINT_THROW("enqueueing edge for addition", err);
            } else
                edge_hash.erase(e.end[0], e.end[1]);
        }

        if (verbose)
            cout << edge_q.size() << " put in q" << endl;

        int num_added = 0;

        while ((!finished && edge_q.size() >= PoolSize) ||
               (finished && penalty < -RecoverMaxPen && !edge_q.empty())) {
            std::sort(edge_q.begin(), edge_q.end());

            try {
                vector<Graph::Edge> add_batch = pool_chunk(edge_q);

                num_added = add_batch.size();
                total_added += num_added;
                core_lp.add_edges(add_batch, false);
                core_lp.primal_opt();
            } CMR_CATCH_PRINT_THROW("adding edges to LP and optimizing", err);

            cur_stat = core_lp.get_stat();

            if (verbose)
                cout << "Added " << num_added << ", total "
                     << total_added
                     << " edges to LP, optimized with stat "
                     << cur_stat << endl;

            if (cur_stat != SolStat::Infeas)
                return true;

            try {
                reg_duals.reset();
                price_edges(edge_q, reg_duals, false);
            } CMR_CATCH_PRINT_THROW("getting new duals/re-pricing", err);

            penalty = 0.0;

            edge_q.erase(std::remove_if(edge_q.begin(), edge_q.end(),
                                        [&penalty](const d_PrEdge &e)
                                        -> bool {
                                            if (e.redcost < 0.0)
                                                penalty += e.redcost;
                                            return e.redcost > - Eps::Zero;
                                        }),
                         edge_q.end());

            if (verbose)
                cout << "Reset duals, q now has size "
                     << edge_q.size() << endl;

            if (num_added > 0) {
                start = end1;
                end2 = start + 1;

                finished = false;
                edge_hash.clear();
                try {
                    for (const d_PrEdge &e : edge_q)
                        edge_hash.add(e.end[0], e.end[1], 1);
                } CMR_CATCH_PRINT_THROW("adding q edges to hash", err);
            }
        }
    }

    if (verbose)
        cout << "Added " << total_added << " edges to LP, still infeasible"
             << endl;
    return false;
}

bool Pricer::feas_gen_edges(std::vector<d_PrEdge> &price_elist,
                            int start_ind, int &end1, int &end2)
{
    runtime_error err("Problem in Pricer::feas_gen_edges");

    if (!reg_duals)
        throw runtime_error("Called Pricer::feas_gen_edges without reg duals");

    vector<double> &npi_est = reg_duals->node_pi_est;

    int ncount = core_graph.node_count();
    int i = end1;
    int j = end2;

    if (i >= ncount) {
        i = 0;
        j = 1;
    }

    for (; j < ncount; ++j) {
        if (npi_est[i] + npi_est[j] > 0.0) {
            int hashval = edge_hash.get_val(i, j);
            if (hashval == -1 && core_graph.find_edge_ind(i, j) == -1) {
                try {
                    price_elist.emplace_back(i, j);
                    edge_hash.add(i, j, 1);
                } CMR_CATCH_PRINT_THROW("adding new cand edge", err);
            }

            if (price_elist.size() == gen_max) {
                end1 = i;
                end2 = j + 1;
                return false;
            }
        }
    }

    while ((i = (i + 1) % ncount) != start_ind) {
        for (j = i + 1; j < ncount; j++) {
            if (npi_est[i] + npi_est[j] > 0.0) {
                int hashval = edge_hash.get_val(i, j);
                if (hashval == -1 && core_graph.find_edge_ind(i, j) == -1) {
                    try {
                        price_elist.emplace_back(i, j);
                        edge_hash.add(i, j, 1);
                    } CMR_CATCH_PRINT_THROW("adding new cand edge", err);
                }

                if (price_elist.size() == gen_max) {
                    end1 = i;
                    end2 = j + 1;
                    return false;
                }
            }
        }
    }

    return true;
}

f64 Pricer::exact_lb(bool full)
{
    vector<PrEdge<f64>> priced_edges;
    return exact_lb(full, priced_edges);
}

f64 Pricer::exact_lb(bool full,
                     vector<PrEdge<f64>> &priced_edges)
{
    for (const Sep::HyperGraph &H : ext_cuts.get_cuts())
        if (H.cut_type() == CutType::Non)
            throw logic_error("Pricer::exact_lb was called w non cut present");

    runtime_error err("Problem in Pricer::exact_lb");

    try {
        ex_duals = util::make_unique<LP::DualGroup<f64>>(true, core_lp,
                                                         ext_cuts);
    } CMR_CATCH_PRINT_THROW("constructing exact DualGroup", err);

    priced_edges.clear();

    vector<double> d_pi;
    vector<char> senses;
    vector<double> rhs_vec;

    int numrows = core_lp.num_rows();
    int numcols = core_lp.num_cols();
    int ncount = inst.node_count();

    vector<f64> ex_pi;

    try {
        core_lp.get_pi(d_pi, 0, numrows - 1);
        senses = core_lp.senses(0, numrows - 1);
        core_lp.get_rhs(rhs_vec, 0, numrows -1);

        ex_pi = vector<f64>(d_pi.begin(), d_pi.end());
    } CMR_CATCH_PRINT_THROW("grabbing cut info", err);

    for (int i = ncount; i < numrows; ++i)
        if (senses[i] == 'G') {
            if (ex_pi[i] < 0.0)
                ex_pi[i] = 0.0;
            else if (senses[i] == 'L')
                if (ex_pi[i] < 0.0)
                    ex_pi[i] = 0.0;
        }

    f64 pi_sum{0.0};

    for (int i = 0; i < numrows; ++i)
        util::add_mult(pi_sum, ex_pi[i], rhs_vec[i]);


    vector<PrEdge<f64>> target_edges;

    if (full) {
        for (int i = 0; i < ncount; ++i)
            for (int j = 0; j < ncount; ++j)
                target_edges.emplace_back(i, j);
    } else {
        for (const Graph::Edge &e : core_graph.get_edges())
            target_edges.emplace_back(e.end[0], e.end[1]);
    }

    price_edges(target_edges, ex_duals, true);

    f64 rc_sum{0.0};

    for (const PrEdge<f64> &e : target_edges)
        if (e.redcost < 0.0)
            rc_sum -= e.redcost;

    f64 bound = pi_sum - rc_sum;

    priced_edges = std::move(target_edges);

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

    // return bound;
}

void Pricer::elim_edges(bool make_opt)
{
    runtime_error err("Problem in Pricer::elim_edges.");

    if (make_opt) {
        double ot = util::zeit();
        cout << "Optimizing before elimination...";
        try { core_lp.primal_opt(); } CMR_CATCH_PRINT_THROW("optimizing", err);
        cout << "obj val " << core_lp.get_objval() << " in "
             << (util::zeit() - ot) << "s" << endl;
    }

    f64 tourlen{core_lp.global_ub()};
    f64 lower_bd{0.0};
    vector<PrEdge<f64>> graph_edges;

    try {
        lower_bd = exact_lb(false, graph_edges);
    } CMR_CATCH_PRINT_THROW("getting exact lb", err);

    f64 gap{tourlen - lower_bd};
    f64 cutoff{gap - 1};

    if (cutoff < 0) {
        if (verbose)
            cout << "Negative cutoff, do not elim." << endl;
        return;
    }

    if (verbose)
        cout << "\tElimination cutoff " << cutoff << endl;

    try {
        util::ptr_reset(ex_duals, true, core_lp, core_lp.external_cuts());
    } CMR_CATCH_PRINT_THROW("getting exact duals", err);

    vector<int> col_delset;
    int ecount = core_graph.edge_count();

    try {
        col_delset.resize(ecount, 0);
    } CMR_CATCH_PRINT_THROW("reserving/prepping vectors", err);

    const vector<int> &tour_colstat = core_lp.get_active_tour().base().colstat;
    const vector<double> &tour_edges = core_lp.get_active_tour().edges();

    int elimct = 0;

    for (int i = 0; i < ecount; ++i) {
        const PrEdge<f64> &e = graph_edges[i];
        if (tour_edges[i] != 0.0 || tour_colstat[i] != 0)
            continue;

        if (e.redcost > cutoff) {
            col_delset[i] = 1;
            ++elimct;
        }
    }

    cout << "\t\t" << elimct << " edges can be eliminated\n" << endl;
    if (elimct == 0)
        return;

    try {
        core_lp.remove_edges(std::move(col_delset));
    } CMR_CATCH_PRINT_THROW("adjusting CoreLP/CoreGraph stuff", err);
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
    const Graph::AdjList &alist = core_graph.get_adj();

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
