/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Interface to the LP solver.
 * The class and prototypes in this file define operations that need to be
 * supported by an LP solver for use with Camargue, with an opaque
 * implementation pointer for the solver being used. Thus, this header is
 * open for re-implementation in other solvers, although only CPLEX is
 * supported for now.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_LP_INTERFACE_H
#define CMR_LP_INTERFACE_H

#include "lp_util.hpp"

#include <utility>
#include <vector>
#include <memory>

namespace CMR {

namespace Sep {
class MIRgroup;
}

namespace LP {

/** Class for storing an lp relaxation via interface to an lp solver.
 * This structure stores the environment and lp data structures of the lp
 * solver, and contains methods for modifying the relaxation in the sense of
 * adding/removing constraints and variables.
 */
class Relaxation {
public:
    /**@name Constructors, destructor, and assignment operators. */
    ///@{


    Relaxation(); //!< Default construct with empty, valid resource handles.

    Relaxation(Relaxation &&lp) noexcept; //!< Move construct.
    Relaxation& operator=(Relaxation &&lp) noexcept;  //!< Move assign.

    Relaxation(const Relaxation &lp) = delete;
    Relaxation& operator=(const Relaxation &lp) = delete;

    ~Relaxation(); //!< Destruct and free resource handles.

    ///@}

    /**@name Methods for querying the relaxation. */
    ///@{

    bool primal_feas() const; //!< Is the resident basis primal feasible.
    bool dual_feas() const; //!< Is the resident basis dual feasible.

    int num_rows() const; //!< Number of rows in the model.
    int num_cols() const; //!< Number of columsn in the model.

    int it_count() const; //!< Number of iterations from last optimization.

    double get_coeff(int row, int col) const; //!< Get constraint matrix entry.

    /// Get a ranged vector of righthand sides.
    void get_rhs(std::vector<double> &rhs, int begin, int end) const;

    /// Get the lower bounds for a range of variables.
    std::vector<double> lower_bds(int begin, int end) const;

    /// Get the upper bounds for a range of variables.
    std::vector<double> upper_bds(int begin, int end) const;

    /// Get the senses of a range of constraints.
    std::vector<char> senses(int begin, int end) const;

    /// Get a sparse coefficient column.
    void get_col(int col, std::vector<int> &cmatind,
                 std::vector<double> &cmatval) const;

    /// Get a sparse coefficient row.
    void get_row(int row, std::vector<int> &rmatind,
                 std::vector<double> &rmatval) const;

    SparseRow get_row(int row) const; //!< Return a sparse coefficient row.

    /// Feasibility status of a given solution with respect to rows.
    void get_row_infeas(const std::vector<double> &x,
                        std::vector<double> &feas_stat, int begin,
                        int end) const;

    /// Feasibility status of a given solution with respect to column bounds.
    void get_col_infeas(const std::vector<double> &x,
                        std::vector<double> &feas_stat, int begin,
                        int end) const;

    double get_objval() const; //!< Objective value for resident solution.

    SolStat get_stat() const; //!< Get the SolStat from last optimization.

    double condition_num() const; //!< Condition number for resident solution.

    void get_x(std::vector<double> &x) const; //!< Get current solution.

    std::vector<double> lp_vec() const; //!< Return current solution.


    /// Get the resident basis.
    void get_base(std::vector<int> &colstat, std::vector<int> &rowstat) const;

    Basis basis_obj() const; //!< Return the resident basis.

    std::vector<int> col_stat() const; //!< Return the column basis.
    std::vector<int> row_stat() const; //!< Return the row basis.

    /// Get constraint slacks for resident solution.
    void get_row_slacks(std::vector<double> &slack, int begin,
                        int end) const;

    /// Return constraint slacks for the resident solution.
    std::vector<double> row_slacks(int begin, int end) const;

    /// Get a range of dual values.
    void get_pi(std::vector<double> &pi, int begin, int end) const;

    /// Return a range of dual values.
    std::vector<double> pi(int begin, int end) const;

    /// Get a range of reduced costs.
    void get_redcosts(std::vector<double> &redcosts, int begin, int end) const;

    /// Return a range of reduced costs.
    std::vector<double> redcosts(int begin, int end) const;

    ///@}

    /**@name Branching methods. */
    ///@{

    /// Get strong branching style estimates for variables.
    void primal_strong_branch(const std::vector<double> &tour_vec,
                              const std::vector<int> &colstat,
                              const std::vector<int> &rowstat,
                              const std::vector<int> &indices,
                              std::vector<Estimate> &down_est,
                              std::vector<Estimate> &up_est,
                              std::vector<Basis> &contra_bases,
                              int itlim, double upperbound);

    /// Tighten the bound on a variable.
    void tighten_bound(int index, char sense, double val);

    /// Change an objective function coefficient.
    void change_obj(int index, double val);

    ///@}


    /**@name Row and column operations. */
    ///@{


    void new_row(char sense, double rhs); //!< Create a new empty row.

    void new_rows(const std::vector<char> &sense,
                  const std::vector<double> &rhs); //!< Create new empty rows.

    void add_cut(double rhs, char sense,
                 const std::vector<int> &rmatind,
                 const std::vector<double> &rmatval); //!< Add constraint row.

    void add_cut(const SparseRow &sp_row); //!< Add constraint row struct.

    void add_cuts(const std::vector<double> &rhs,
                  const std::vector<char> &sense,
                  const std::vector<int> &rmatbeg,
                  const std::vector<int> &rmatind,
                  const std::vector<double> &rmatval); //!< Add constraint rows.

    /// Delete a not-necessarily-contiguous set of rows.
    void del_set_rows(std::vector<int> &delstat);

    void add_col(const double objval,
                 const std::vector<int> &indices,
                 const std::vector<double> &coeffs,
                 const double lb, const double ub); //!< Add a column.

    /// Delete a not-necessarily-contiguous set of columns.
    void del_set_cols(std::vector<int> &delstat);

    ///@}

    /**@name Basis operations. */
    ///@{

    /// Instate a basis in the problem corresponding to a solution.
    void copy_start(const std::vector<double> &x);

    /// Instate a specified basis in the problem, with specified solution.
    void copy_start(const std::vector<double> &x,
                    const std::vector<int> &col_stat,
                    const std::vector<int> &row_stat);

    /// Instate a basis without explicitly specifying a solution.
    void copy_base(const std::vector<int> &col_stat,
                   const std::vector<int> &row_stat);

    /// Overload of copy_base for basis struct.
    void copy_base(const Basis &base);

    void factor_basis(); //!< Make the current basis resident in the problem.

    ///@}

    /**@name Pivoting and optimizing. */
    ///@{

    void switch_steepest(); //!< Permanent switch to steepest edge pivots.

    void primal_opt(); //!< Optimize the Relaxation with primal simplex.
    void dual_opt(); //!< Optimize the relaxation with dual simplex.
    void nondegen_pivot(double upper_bound); //!< Primal non-degenerate pivot.

    /// Do a primal non-degenerate pivot with a callback to record a basis.
    void cb_nondegen_pivot(double upper_bound, Basis &base, int bas_freq);

    void one_primal_pivot(); //!< Perform exactly one primal simplex pivot.
    void one_dual_pivot(); //!< Perform exactly one dual simplex pivot.

    void primal_recover(); //!< Pivot until the basis is primal feasible.

    ///@}

    void init_mir_data(Sep::MIRgroup &mir_data); //!< Construct a Sep::MirGroup.

private:
    struct solver_impl; //!< Implementation hiding for solver.
    std::unique_ptr<solver_impl> simpl_p; //!< Pointer to solver implementation.
};

/// Handle for data during a non-degenerate pivot callback.
struct NDpivotHandle {
    NDpivotHandle(Relaxation &_rel, double ub, int bfreq)
        : rel(_rel), tour_base(), upper_bound(ub), basis_freq(bfreq) {}

    Relaxation &rel;
    Basis tour_base; //!< A new basis for the tour to be set by the callback.

    double upper_bound; //!< Value of best primal feasible solution.

    /// How often to copy a basis.
    /// Zero for only the first primal feasible basis.
    /// One for every new primal feasible basis.
    /// `n > 1` for only when `pfeas_itcount % n == 0`
    int basis_freq;

    int pfeas_itcount = 0; //!< Number of pivots since becoming primal feasible.
};

}
}

#endif
