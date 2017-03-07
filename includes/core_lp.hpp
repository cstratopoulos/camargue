/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Managing Core LP relaxations of TSP instances.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CORE_LP_H
#define CMR_CORE_LP_H

#include "lp_interface.hpp"
#include "active_tour.hpp"
#include "cc_lpcuts.hpp"
#include "process_cuts.hpp"
#include "datagroups.hpp"
#include "hypergraph.hpp"
#include "cutmon.hpp"
#include "err_util.hpp"
#include "util.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CMR {

class Solver;

namespace LP {


/** Class for storing the core lp associated to a TSP instance and pivoting.
 * This class contains the edges and constraints currently under consideration
 * in an lp relaxation of a TSP instance.
 */
class CoreLP : public Relaxation {
public:
    CoreLP(Graph::CoreGraph &core_graph_,
           Data::BestGroup &best_data_);

    /// Compute a primal non-degenerate pivot from the active_tour.
    LP::PivType primal_pivot();

    void pivot_back(bool prune_slacks); //!< Pivot back to active_tour.

    void add_cuts(Sep::LPcutList &cutq);
    void add_cuts(Sep::CutQueue<Sep::dominoparity> &dp_q);
    void add_cuts(Sep::CutQueue<LP::SparseRow> &gmi_q);
    void add_cuts(Sep::CutQueue<Sep::ex_blossom> &ex2m_q);
    void add_cuts(Sep::CutQueue<Sep::HyperGraph> &pool_q);

    void add_edges(const std::vector<Graph::Edge> &add_batch,
                   bool reinstate);

    void remove_edges(std::vector<int> edge_delstat);


    /// Returns true iff x_vec satisfies all constraints and column bounds.
    template<typename numtype>
    bool check_feas(const std::vector<numtype> &x_vec);

    /// Get a const reference to the SupportGroup for the most recent pivot.
    const Data::SupportGroup &support_data() const { return supp_data; }

    const Sep::ExternalCuts &external_cuts() const { return ext_cuts; }

    const LP::CutMonitor &cut_monitor() const { return cut_mon; }

    const LP::ActiveTour &get_active_tour() const { return active_tour; }

    void set_active_tour(std::vector<int> tour_nodes);

    /// Average number of iterations per primal_pivot.
    int avg_itcount() const { return sum_it_count / num_nd_pivots; }

    double active_tourlen() const { return active_tour.length(); }

    double global_ub() const { return best_data.min_tour_value; }


    friend class CMR::Solver;

private:
    /// Instate active_tour, checking feasibility if it throws.
    void instate_active();

    /// As above but with a reset instate.
    void reset_instate_active();

    void handle_aug_pivot(std::vector<int> tour_nodes, Basis aug_base);

    void prune_slacks();

    void rebuild_basis();

    void purge_gmi();

    Graph::CoreGraph &core_graph;
    Data::BestGroup &best_data;
    Data::SupportGroup supp_data;

    Sep::ExternalCuts ext_cuts;

    LP::CutMonitor cut_mon;

    ActiveTour active_tour;

    std::vector<double> lp_edges;
    std::vector<double> pi_vals;
    std::vector<double> feas_stat;

    int num_nd_pivots = 0;
    int sum_it_count = 0;

    int prev_numrows;
};

//////////////////// TEMPLATE METHOD IMPLEMENTATIONS //////////////////////////

template<typename numtype>
bool CoreLP::check_feas(const std::vector<numtype> &x_vec)
{
    using std::cout;
    using std::endl;
    using std::runtime_error;
    using std::vector;
    namespace Eps = Epsilon;

    runtime_error err("Problem in CoreLP::check_feas");
    bool result = true;

    int numrows = num_rows();
    int numcols = num_cols();

    vector<double> dbl_x;
    vector<double> row_feas;
    vector<double> col_feas;

    cout << "Reporting if solution is feasible, " << numrows << " rows and "
         << numcols << " cols in the LP" << endl;

    try {
        dbl_x = vector<double>(x_vec.begin(), x_vec.end());
        get_row_infeas(dbl_x, row_feas, 0, numrows - 1);
        get_col_infeas(dbl_x, col_feas, 0, numcols - 1);
    } CMR_CATCH_PRINT_THROW("making solver queries", err);

    int ncount = core_graph.node_count();

    for (int i = 0; i < numrows; ++i) {
        double rowfeas = row_feas[i];
        if (rowfeas != 0.0) {
            result = false;
            cout << "\nRow " << i << " has infeas of "  << rowfeas << " on ";
            if (i < ncount)
                cout << "degree equation" << endl;
            else {
                const Sep::HyperGraph &H = ext_cuts.get_cut(i);
                cout << H.cut_type() << ", sense "
                     << H.get_sense() << ", rhs " << H.get_rhs() << endl;
                LP::SparseRow rel_row;
                LP::SparseRow hg_row;

                try {
                    rel_row = get_row(i);
                    H.get_coeffs(core_graph.get_edges(), hg_row.rmatind,
                                 hg_row.rmatval);
                } CMR_CATCH_PRINT_THROW("getting SparseRows", err);
                cout << "HG row act: " << Sep::get_activity(dbl_x, hg_row)
                     << ", rel row act: " << Sep::get_activity(dbl_x, rel_row)
                     << endl;
                cout << "HG row size " << hg_row.rmatind.size() << ", rel "
                     << rel_row.rmatind.size() << endl;
                cout << "Rel row sense " << rel_row.sense << ", rhs "
                     << rel_row.rhs << endl;

                cout << "Checking HG row for zeros or fracs...." << endl;
                for (int i = 0; i < hg_row.rmatind.size(); ++i) {
                    int ind = hg_row.rmatind[i];
                    double val = hg_row.rmatval[i];

                    if (val == 0 || val != static_cast<int>(val))
                        cout << "Index " << ind << " has val " << val << "\n";

                }

                cout << "Checking rel row for zeros or fracs...." << endl;
                for (int i = 0; i < rel_row.rmatind.size(); ++i) {
                    int ind = rel_row.rmatind[i];
                    double val = rel_row.rmatval[i];

                    if (val == 0 || val != static_cast<int>(val))
                        cout << "Index " << ind << " has val " << val << "\n";

                }
            }
        }
    }

    for (int i = 0; i < numcols; ++i) {
        double colfeas = col_feas[i];
        if (colfeas != 0.0) {
            result = false;
            cout << "Found infeas of " << colfeas << " on edge "
                 << core_graph.get_edge(i) << endl;
        }
    }

    cout << "Report completed with result " << result << endl;
    return result;
}


}
}

#endif
