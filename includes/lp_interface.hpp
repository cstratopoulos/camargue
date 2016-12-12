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

/** Class for storing the lp relaxation of a TSP instance. 
 * This structure stores the environment and lp data structures of the lp 
 * solver, and contains methods for modifying the relaxation in the sense of
 * adding/removing constraints and variables.
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

    /**@name Row and column queries/operations. */
    ///@{

    int num_rows();
    int num_cols();

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

    void get_row_infeas(const std::vector<double> &x,
                        std::vector<double> &feas_stat, int begin, int end);

    void del_set_rows(std::vector<int> &delstat);

    void new_cols(const std::vector<double> &objfun,
                  const std::vector<int> &cmatbeg,
                  const std::vector<int> &cmatind,
                  const std::vector<double> &cmatval,
                  const std::vector<double> &lb,
                  const std::vector<double> &ub);
    
    
    ///@}

    /**@name Basis queries/operations. */
    ///@{

    void get_base(std::vector<int> &colstat, std::vector<int> &rowstat);
    
    void copy_start(const std::vector<double> &x);
    void copy_start(const std::vector<double> &x,
                    const std::vector<int> &col_stat,
                    const std::vector<int> &row_stat);

    void factor_basis();
    
    ///@}

    /**@name Pivoting and optimizing. */
    ///@{

    void nondegen_pivot(const double lowlimit);
    
    ///@}

    /**@name Solution queries. */
    ///@{

    double get_objval();
    void get_x(std::vector<double> &x);

    bool dual_feas();

    void get_row_slacks(std::vector<double> &slack, int begin, int end);
    void get_pi(std::vector<double> &pi, int begin, int end);
    
    ///@}

private:
    CPXENVptr cplex_env;
    CPXLPptr cplex_lp;
};

}
}

#endif
