/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief LP SOLVER INTERFACE STRUCT AND METHODS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_LP_INTERFACE_H
#define CMR_LP_INTERFACE_H

#include <vector>
#include <memory>

namespace CMR {
namespace LP {

enum BStat {
    AtLower = 0,
    Basic = 1,
    AtUpper = 2,
    FreeSuper = 3
};

struct Basis {
    Basis() = default;
    
    std::vector<int> colstat;
    std::vector<int> rowstat;
};


/** Class for storing an lp relaxation via interface to an lp solver. 
 * This structure stores the environment and lp data structures of the lp 
 * solver, and contains methods for modifying the relaxation in the sense of
 * adding/removing constraints and variables.
 */
class Relaxation {
public:
    /**@name Constructors, destructor, and assignment operators. */
    ///@{

    
    Relaxation(); /**< Default constructor. */
    
    Relaxation(Relaxation &&lp) noexcept; /**< Move constructor. */
    Relaxation& operator=(Relaxation &&lp) noexcept; /**< Move assignment. */

    Relaxation(const Relaxation &lp) = delete;
    Relaxation& operator=(const Relaxation &lp) = delete;
    
    ~Relaxation(); /**< Destruct and free associated memory. */
    
    ///@}

    /**@name Methods for querying the relaxation. */
    ///@{

    bool dual_feas() const; /**< Is the resident basis dual feasible. */

    int num_rows() const; /**< Number of rows in the model. */
    int num_cols() const; /**< Number of columsn in the model. */

    int it_count() const;

    /** Check the feasibility status of a given solution.
     * If \p x is an lp solution, the `i`th entry of \p feas_stat will be
     * nonzero if \p x violates the constraint in row `i`.
     */
    void get_row_infeas(const std::vector<double> &x,
                        std::vector<double> &feas_stat, int begin,
                        int end) const;

    double get_objval() const; /**< Objective value for resident solution. */
    
    void get_x(std::vector<double> &x) const; /**< Get current solution. */
    
    std::vector<double> lp_vec() const; /**< Return current solution. */

    /** Get the resident basis. */
    void get_base(std::vector<int> &colstat,
                  std::vector<int> &rowstat) const;

    /** Get constraint slacks for the resident solution. */
    void get_row_slacks(std::vector<double> &slack, int begin,
                        int end) const;

    /** Return constraint slacks for the resident solution. */
    std::vector<double> row_slacks(int begin, int end) const;

    /** Get a range of dual values. */
    void get_pi(std::vector<double> &pi, int begin, int end) const;

    /** Return a range of dual values. */
    std::vector<double> pi(int begin, int end) const;

    /// Get a range of reduced costs.
    void get_redcosts(std::vector<double> &redcosts, int begin, int end) const;

    /// Return a range of reduced costs.
    std::vector<double> redcosts(int begin, int end) const;

    ///@}

    /**@name Branch variable selection methods. */
    ///@{

    /** Get the Driebeek penalties for specified indices. */
    void get_penalties(const std::vector<int> &indices,
                       std::vector<double> &downratio,
                       std::vector<double> &upratio);

    /** Get strong branching objective values. 
     * This function will evaluate the variables in \p indices (which should
     * be the indices of fractional basic variables) as candidates for 
     * branching edges by performing at most \p itlim dual steepest edge
     * simplex pivots, with an objective function upper bound of
     * \p upperbound.  If \p opt_first is true, the function will optimize
     * the relaxation first, saving the resident basis/solution which will
     * be restored afterwards. Otherwise, the pivots will be performed 
     * from the resident basis/solution. 
     */
    void dual_strong_branch(const std::vector<int> &indices,
                            std::vector<double> &downobj,
                            std::vector<double> &upobj, int itlim,
                            double upperbound);

    /** Get modified strong branch objective values.
     * This function behaves like dual_strong_branch, but with
     * an implementation particular to this project. Here \p tour_vec shall
     * be a binary vector indicating edges in a tsp tour, with associated
     * column and row basis provided by \p colstat and \p rowstat, 
     * respectively. Instead of using dual
     * steepest edge simplex pivots, this function will perform at most
     * \p itlim primal steepeste edge simplex pivots, with the tour solution
     * and  basis as a starting point. 
     */
    void primal_strong_branch(const std::vector<double> &tour_vec,
                              const std::vector<int> &colstat,
                              const std::vector<int> &rowstat,
                              const std::vector<int> &indices,
                              std::vector<double> &downobj,
                              std::vector<double> &upobj,
                              int itlim, double upperbound);

    ///@}


    /**@name Row and column operations. */
    ///@{

    /** Create a new empty row. */
    void new_row(const char sense, const double rhs);

    /** Create multiple new empty rows. */
    void new_rows(const std::vector<char> &sense,
                  const std::vector<double> &rhs);

    /** Add a constraint row to the model. */
    void add_cut(const double rhs, const char sense,
                 const std::vector<int> &rmatind,
                 const std::vector<double> &rmatval);

    /** Add multiple constraint rows to the model at once. */
    void add_cuts(const std::vector<double> &rhs,
                  const std::vector<char> &sense,
                  const std::vector<int> &rmatbeg,
                  const std::vector<int> &rmatind,
                  const std::vector<double> &rmatval);


    /** Delete a specified, not-necessarily-continuous set of rows. */
    void del_set_rows(std::vector<int> &delstat);

    /** Add a column to the model. */
    void add_col(const double objval,
                 const std::vector<int> &indices,
                 const std::vector<double> &coeffs,
                 const double lb, const double ub);
    
    ///@}

    /**@name Basis operations. */
    ///@{

    /** Instate a basis in the problem corresponding to a solution. */
    void copy_start(const std::vector<double> &x);

    /** Instate a specified basis in the problem, with specified solution. */
    void copy_start(const std::vector<double> &x,
                    const std::vector<int> &col_stat,
                    const std::vector<int> &row_stat);

    /** Make the current basis/solution resident in the problem. */
    void factor_basis();
    
    ///@}

    /**@name Pivoting and optimizing. */
    ///@{

    /** Optimize the relaxation with primal simplex. */
    void primal_opt();

    /** Optimize the relaxation with dual simplex. */
    void dual_opt();

    /** Find a primal non-degenerate pivot. */
    void nondegen_pivot(const double lowlimit);

    /** Perform exactly one primal simplex pivot. */
    void single_pivot();
    
    ///@}

private:
    struct solver_impl;
    std::unique_ptr<solver_impl> simpl_p;
};

}
}

#endif
