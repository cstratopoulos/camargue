/**
 * @file
 * @brief Pricer implementations for recovering an infeasible LP.
 */

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

namespace CMR {

using LP::PivType;
using CutType = Sep::HyperGraph::Type;
using f64 = util::Fixed64;

namespace Eps = Epsilon;

namespace Price {

using d_PrEdge = PrEdge<double>;

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

}
}
