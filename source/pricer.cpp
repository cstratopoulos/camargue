#include "pricer.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <utility>

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

    CCtsp_edgegenerator *current_eg;

    if (piv_stat == PivType::FathomedTour){
        current_eg = &eg_full;
        if (verbose)
            cout << "\tRunning full eg\n";
    } else if (piv_stat == PivType::Tour){
        current_eg = &eg_inside;
        if (verbose)
            cout << "\tRunning inside eg\n";
    } else
        throw logic_error("Tried to run pricing on non tour.");

    if (CCtsp_reset_edgegenerator(current_eg,
                                  &(reg_duals->node_pi_est[0]), !verbose)) {
        cerr << "CCtsp_reset_edgegenerator failed.\n";
        throw err;
    }

    CCrandstate rstate;

    CCutil_sprand(inst.seed(), &rstate);

    int finished = 0;
    int outercount = 0;
    int total_added = 0;

    double penalty = 0.0;

    double upper_bound = core_lp.global_ub();
    //double lower_bound = core_lp.get_objval();

    vector<d_PrEdge> price_elist;
    vector<d_PrEdge> edge_q;

    while (!finished) {
        if (verbose)
            cout << "\tEntering EG loop, pass " << ++outercount << "\n\t";

        int num_gen = 0;

        price_elist.clear();

        if (CCtsp_generate_edges(current_eg, gen_max, &num_gen, &gen_elist[0],
                                 &gen_elen[0], &finished, !verbose, &rstate)) {
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

        if (verbose)
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

        if (verbose)
            cout << "\t" << edge_q.size() << " edges in edge_q\n"
                 << "\t" << penalty << " penalty accrued, finished: "
                 << finished << "\n";

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
        int innercount = 0;

        while (((!finished && edge_q.size() >= PoolSize) ||
                (finished && penalty < -MaxPenalty && !edge_q.empty()))) {
            if (verbose)
                cout << "\t\tOpt tour inner price loop, pass " << ++innercount
                     << "\n";
            std::sort(edge_q.begin(), edge_q.end());

            try {
                vector<Graph::Edge> add_batch = pool_chunk(edge_q);

                num_added = add_batch.size();
                total_added += num_added;
                core_lp.add_edges(add_batch, true);
                core_lp.primal_opt();
                new_objval = core_lp.get_objval();
            } CMR_CATCH_PRINT_THROW("adding edges to lp and optimizing", err);

            if (verbose)
                cout << "\t\tAdded " << num_added
                     << " edges in opt solution.\n";

            if (new_objval <= upper_bound - 1.0 + Eps::Zero) {
                if (innercount == 1)
                    if (verbose)
                        cout << "\tTour no longer optimal after adding edges\n";
                result = ScanStat::Full;
            } else {
                if (verbose)
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
                                              !verbose)) {
                    cerr << "CCtsp_reset_edgegenerator failed.\n";
                    throw err;
                }

                finished = 0;
                edge_hash.clear();
                for (const d_PrEdge &e : edge_q)
                    edge_hash.add(e.end[0], e.end[1], 1);
            }

            if (verbose)
                cout << "\t\tInner penalty: "
                     << penalty << ", new edge_q size "
                     << edge_q.size() << ", finished: " << finished << "\n";
        }
    }

    cout << "added " << total_added << " edges, " << result << "\n\t"
         << "objval " << core_lp.get_objval() << ", dual feas "
         << core_lp.dual_feas() << ", col count "
         << core_lp.num_cols() << endl;

    if (try_elim && result != ScanStat::FullOpt && core_lp.dual_feas()) {
        try { elim_edges(false); }
        CMR_CATCH_PRINT_THROW("eliminating after gen", err);
    }

    return result;
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

    price_edges(target_edges, ex_duals);

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
