/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief LP SOLVER INTERFACE STRUCT AND METHODS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_LP_INTERFACE_H
#define CMR_LP_INTERFACE_H

#include <cplex.h>

#include <vector>

namespace CMR {
namespace LP {

/** Class for storing an lp relaxation via interface to an lp solver. 
 * This structure stores the environment and lp data structures of the lp 
 * solver, and contains methods for modifying the relaxation in the sense of
 * adding/removing constraints and variables.
 * @remark This is implemented as an interface to CPLEX 12, but use of 
 * other solvers could be supported by re-implementing the protected members
 * and changing private members appropriately. 
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

    /** Get the dual values. */
    void get_pi(std::vector<double> &pi, int begin, int end) const;

    /** Return the dual values. */
    std::vector<double> pi(int begin, int end) const;

    /** Get the Driebeek penalties for specified indices. */
    void get_penalties(const std::vector<int> &indices,
                       std::vector<double> &downratio,
                       std::vector<double> &upratio);

    /** Get strong branching objective values. */
    void dual_strong_branch(const std::vector<int> &indices,
                            std::vector<double> &downobj,
                            std::vector<double> &upobj, int itlim);


    ///@}



protected:

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

    /** Find a primal non-degenerate pivot.
     * If there is a solution `x` resident in the problem, and `x` has
     * objective value `lowlimit + CMR::Epsilon::Zero`, this function
     * pivots until reaching a solution which differs from `x`. 
     */
    void nondegen_pivot(const double lowlimit);
    
    ///@}

private:
    CPXENVptr cplex_env;
    CPXLPptr cplex_lp;
};

}
}

#endif
