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
    
    Relaxation();
    Relaxation(Relaxation &&lp) noexcept;
    Relaxation& operator=(Relaxation &&lp) noexcept;

    Relaxation(const Relaxation &lp) = delete;
    Relaxation& operator=(const Relaxation &lp) = delete;
    
    ~Relaxation();
    
    ///@}

    /**@name Methods for querying the relaxation. */
    ///@{

    bool dual_feas() const;


    int num_rows() const;
    int num_cols() const;

    void get_row_infeas(const std::vector<double> &x,
                        std::vector<double> &feas_stat, int begin,
                        int end) const;

    double get_objval() const;
    
    void get_x(std::vector<double> &x) const;
    std::vector<double> lp_vec() const;

    void get_base(std::vector<int> &colstat,
                  std::vector<int> &rowstat) const;


    void get_row_slacks(std::vector<double> &slack, int begin,
                        int end) const;
    std::vector<double> row_slacks(int begin, int end) const;
    
    void get_pi(std::vector<double> &pi, int begin, int end) const;
    std::vector<double> pi(int begin, int end) const;

    void get_penalties(const std::vector<int> &indices,
                       std::vector<double> &downratio,
                       std::vector<double> &upratio);


    ///@}



protected:

    /**@name Row and column operations. */
    ///@{

    void new_row(const char sense, const double rhs);
    void new_rows(const std::vector<char> &sense,
                  const std::vector<double> &rhs);
    
    void add_cut(const double rhs, const char sense,
                 const std::vector<int> &rmatind,
                 const std::vector<double> &rmatval);
    
    void add_cuts(const std::vector<double> &rhs,
                  const std::vector<char> &sense,
                  const std::vector<int> &rmatbeg,
                  const std::vector<int> &rmatind,
                  const std::vector<double> &rmatval);


    void del_set_rows(std::vector<int> &delstat);

    void add_col(const double objval,
                 const std::vector<int> &indices,
                 const std::vector<double> &coeffs,
                 const double lb, const double ub);
    
    ///@}

    /**@name Basis operations. */
    ///@{

    
    void copy_start(const std::vector<double> &x);
    void copy_start(const std::vector<double> &x,
                    const std::vector<int> &col_stat,
                    const std::vector<int> &row_stat);

    void factor_basis();
    
    ///@}

    /**@name Pivoting and optimizing. */
    ///@{

    void primal_opt();

    void nondegen_pivot(const double lowlimit);
    
    ///@}

private:
    CPXENVptr cplex_env;
    CPXLPptr cplex_lp;
};

}
}

#endif
