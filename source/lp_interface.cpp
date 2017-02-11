#include "lp_interface.hpp"
#include "err_util.hpp"
#include "util.hpp"

#if CMR_HAVE_SAFEGMI

#include <safemir/src/cplex_slvr.cpp>
#include <safemir/src/ds_slvr.cpp>

#endif

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <string>
#include <typeinfo>
#include <utility>


#include <cplex.h>

using std::abs;

using std::vector;

using std::unique_ptr;

using std::cout;
using std::cerr;
using std::endl;
using std::to_string;
using std::string;

using std::runtime_error;
using std::logic_error;
using std::exception;


using cpx_err = CMR::util::retcode_error;


namespace CMR {

namespace Eps = CMR::Epsilon;

namespace LP {

constexpr double CPXzero = 1E-10;
constexpr double CPXint_tol = 0.0001;

struct Relaxation::solver_impl {
    solver_impl();
    ~solver_impl();
    
    CPXENVptr env;
    CPXLPptr lp;
};

Relaxation::solver_impl::solver_impl() try
{
    int rval = 0;
    
    lp = (CPXLPptr) NULL;
    env = CPXopenCPLEX(&rval);

    if (rval) 
        throw cpx_err(rval, "CPXopenCPLEX");

    // this scope guard should prevent memory leaks in the event of an error
    // before the actual problem is initialized. 
    auto cleanup = util::make_guard([&rval, this] {
        if (rval)
            CPXcloseCPLEX(&env);
    });

    rval = CPXsetintparam(env, CPX_PARAM_THREADS, 1);
    if (rval)
        throw cpx_err(rval, "CPXsetintparam thread count");

    rval = CPXsetintparam(env, CPX_PARAM_PERIND, CPX_OFF);
    if (rval)
        throw cpx_err(rval, "CPXsetintparam perturbation");

    rval = CPXsetintparam(env, CPX_PARAM_AGGIND, CPX_OFF);
    if (rval)
        throw cpx_err(rval, "CPXsetintparam aggregator");
    
    rval = CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
    if (rval)
        throw cpx_err(rval, "CPXsetintparam presolve");

    rval = CPXsetintparam(env, CPX_PARAM_PPRIIND, CPX_PPRIIND_DEVEX);
    if (rval) 
        throw cpx_err(rval, "CPXsetintparam primal pricing");

    rval = CPXsetintparam(env, CPX_PARAM_DPRIIND, CPX_DPRIIND_STEEP);
    if (rval)
        throw cpx_err(rval, "CPXsetintparam dual pricing");
    
    char unused;

    lp = CPXcreateprob(env, &rval, &unused);

    if (rval) {
        throw cpx_err(rval, "CPXcreateprob");
    }
} catch (const exception &e) {
    throw runtime_error("cplex solver_impl constructor failed.");
}

Relaxation::solver_impl::~solver_impl()
{
    if (env) {
        if (lp) {
            CPXfreeprob(env, &lp);
            lp = (CPXLPptr) NULL;
        }
        CPXcloseCPLEX(&env);
        env = (CPXENVptr) NULL;
    }
}


template <typename cplex_query>
std::vector<double> info_vec(cplex_query F, const char* Fname,
                             const CPXENVptr cplex_env,
                             const CPXLPptr cplex_lp,
                             int begin, int end)
{
    std::vector<double> result(end - begin + 1);

    int rval = F(const_cast<CPXENVptr>(cplex_env),
                 const_cast<CPXLPptr>(cplex_lp), &result[0], begin, end);
    if (rval)
        throw cpx_err(rval, Fname);

    return result;
}

template <typename cplex_query, typename vectype>
void set_info_vec(cplex_query F, const char *Fname,
                  const CPXENVptr cplex_env, const CPXLPptr cplex_lp,
                  vector<vectype> &info_vec, int begin, int end)
{
    info_vec.resize(end - begin + 1);

    int rval = F(cplex_env, cplex_lp, &info_vec[0], begin, end);
    if (rval)
        throw cpx_err(rval, Fname);
}

static int pfeas_cb (CPXCENVptr cpx_env, void *cbdata, int wherefrom,
                     void *cbhandle)
{
    int pfeas = 0;
    
    int rval = CPXgetcallbackinfo(cpx_env, cbdata, wherefrom,
                                  CPX_CALLBACK_INFO_PRIMAL_FEAS,
                                  &pfeas);
    if (rval)
        throw cpx_err(rval, "CPXgetcallbackinfo pfeas.");
    
    return pfeas;
}

class CPXintParamGuard {
public:
    CPXintParamGuard(int which, int new_value, CPXENVptr env,
                  const string &p_desc) try
        : which_param(which), cplex_env(env), param_desc(p_desc)
    {
        int rval = CPXgetintparam(cplex_env, which_param, &old_value);
        
        if (rval)
            throw cpx_err(rval, "CPXgetintparam " + param_desc);

        rval = CPXsetintparam(env, which_param, new_value);
        if (rval)
            throw cpx_err(rval, "CPXsetintparam " + param_desc);
    } catch (const exception & e) {
        cerr << e.what() << "\n";
        throw runtime_error("Couldn't create int param guard.");
    }

    ~CPXintParamGuard()
    {
        int rval = CPXsetintparam(cplex_env, which_param, old_value);
        if (rval) {
            cerr << "\tFATAL ERROR: Failed to revert int param "
                 << param_desc << ", rval " << rval << "in destructor.\n";
            exit(1);
        }
    }

private:
    const int which_param;
    CPXENVptr cplex_env;
    const string &param_desc;
    int old_value;
};

class CPXdblParamGuard {
public:
    CPXdblParamGuard(int which, double new_value, CPXENVptr env,
                  const string &p_desc) try
        : which_param(which), cplex_env(env), param_desc(p_desc)
    {
        int rval = CPXgetdblparam(cplex_env, which_param, &old_value);
        
        if (rval)
            throw cpx_err(rval, "CPXgetdblparam " + param_desc);

        rval = CPXsetdblparam(env, which_param, new_value);
        if (rval)
            throw cpx_err(rval, "CPXsetdblparam " + param_desc);
    } catch (const exception & e) {
        cerr << e.what() << "\n";
        throw runtime_error("Couldn't create double param guard.");
    }

    ~CPXdblParamGuard()
    {
        int rval = CPXsetdblparam(cplex_env, which_param, old_value);
        if (rval) {
            cerr << "\tFATAL ERROR: Failed to revert double param "
                 << param_desc << ", rval " << rval << "in destructor.\n";
            exit(1);
        }
    }

private:
    const int which_param;
    CPXENVptr cplex_env;
    const string &param_desc;
    double old_value;
};

class CPXlongParamGuard {
public:
    CPXlongParamGuard(int which, CPXLONG new_value, CPXENVptr env,
                  const string &p_desc) try
        : which_param(which), cplex_env(env), param_desc(p_desc)
    {
        int rval = CPXgetlongparam(cplex_env, which_param, &old_value);
        
        if (rval)
            throw cpx_err(rval, "CPXgetlongparam " + param_desc);

        rval = CPXsetlongparam(env, which_param, new_value);
        if (rval)
            throw cpx_err(rval, "CPXsetlongparam " + param_desc);
    } catch (const exception & e) {
        cerr << e.what() << "\n";
        throw runtime_error("Couldn't create long param guard.");
    }

    ~CPXlongParamGuard()
    {
        int rval = CPXsetlongparam(cplex_env, which_param, old_value);
        if (rval) {
            cerr << "\tFATAL ERROR: Failed to revert long param "
                 << param_desc << ", rval " << rval << "in destructor.\n";
            exit(1);
        }
    }

private:
    const int which_param;
    CPXENVptr cplex_env;
    const string &param_desc;
    CPXLONG old_value;
};


Relaxation::Relaxation()
try : simpl_p(util::make_unique<solver_impl>())
{} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Relaxation constructor failed.");
}

Relaxation::Relaxation(Relaxation &&lp) noexcept
    : simpl_p(std::move(lp.simpl_p))
{
    lp.simpl_p.reset(nullptr);
}

Relaxation& Relaxation::operator=(Relaxation &&lp) noexcept
{
    simpl_p = std::move(lp.simpl_p);
    lp.simpl_p.reset(nullptr);
    
    return *this;
}

Relaxation::~Relaxation() {}

int Relaxation::num_rows() const
{
    return CPXgetnumrows(simpl_p->env, simpl_p->lp);
}

int Relaxation::num_cols() const
{
    return CPXgetnumcols(simpl_p->env, simpl_p->lp);
}

int Relaxation::it_count() const
{
    return CPXgetitcnt(simpl_p->env, simpl_p->lp);
}

double Relaxation::get_coeff(int row, int col) const
{
    double result;
    int rval = CPXgetcoef(simpl_p->env, simpl_p->lp, row, col, &result);
    if (rval)
        throw cpx_err(rval, "CPXgetcoef");
    return result;
}

void Relaxation::get_rhs(vector<double> &rhs, int begin, int end) const
{
    set_info_vec(CPXgetrhs, "CPXgetrhs", simpl_p->env, simpl_p->lp, rhs,
                 begin, end);
}

vector<char> Relaxation::senses(int begin, int end) const
{
    vector<char> result;
    set_info_vec(CPXgetsense, "CPXgetsense", simpl_p->env, simpl_p->lp,
                 result, begin, end);
    return result;
}

void Relaxation::get_col(const int col, vector<int> &cmatind,
                           vector<double> &cmatval) const
{
    int surplus = 0;
    int cmatbeg = 0;
    int nzcnt = 0;

    int rval = CPXgetcols(simpl_p->env, simpl_p->lp,
                          &nzcnt, &cmatbeg,
                          &cmatind[0], &cmatval[0], 0,
                          &surplus, col, col);
    if (rval != CPXERR_NEGATIVE_SURPLUS)
        throw cpx_err(rval, "CPXgetcols initial");

    try {
        cmatind.resize(-surplus);
        cmatval.resize(-surplus);
    } catch (const exception &e) {
        cerr << e.what() << "\n";
        throw runtime_error("Couldn't resize vectors in Relaxation::get_col");
    }

    rval = CPXgetcols(simpl_p->env, simpl_p->lp,
                      &nzcnt, &cmatbeg,
                      &cmatind[0], &cmatval[0], cmatind.size(),
                      &surplus, col, col);
    if (rval)
        throw cpx_err(rval, "CPXgetcols actual");
}

void Relaxation::get_row(const int row, vector<int> &rmatind,
                           vector<double> &rmatval) const
{
    int surplus = 0;
    int rmatbeg = 0;
    int nzcnt = 0;

    int rval = CPXgetrows(simpl_p->env, simpl_p->lp,
                          &nzcnt, &rmatbeg,
                          &rmatind[0], &rmatval[0], 0,
                          &surplus, row, row);
    if (rval != CPXERR_NEGATIVE_SURPLUS)
        throw cpx_err(rval, "CPXgetrows initial");

    try {
        rmatind.resize(-surplus);
        rmatval.resize(-surplus);
    } catch (const exception &e) {
        cerr << e.what() << "\n";
        throw runtime_error("Couldn't resize vectors in Relaxation::get_row");
    }

    rval = CPXgetrows(simpl_p->env, simpl_p->lp,
                      &nzcnt, &rmatbeg,
                      &rmatind[0], &rmatval[0], rmatind.size(),
                      &surplus, row, row);
    if (rval)
        throw cpx_err(rval, "CPXgetrows actual");
}

void Relaxation::new_row(const char sense, const double rhs)
{
    int rval = CPXnewrows(simpl_p->env, simpl_p->lp, 1, &rhs, &sense, NULL,
                          NULL);

    if (rval)
        throw cpx_err(rval, "CPXnewrows");
}

void Relaxation::new_rows(const vector<char> &sense,
                          const vector<double> &rhs)
{
    int rval = CPXnewrows(simpl_p->env, simpl_p->lp, sense.size(), rhs.data(),
                          sense.data(), NULL, NULL);
    if (rval)
        throw cpx_err(rval, "CPXnewrows");
}

void Relaxation::add_cut(const double rhs, const char sense,
                         const vector<int> &rmatind,
                         const vector<double> &rmatval)
{
    int rmatbeg = 0;

    int rval = CPXaddrows(simpl_p->env, simpl_p->lp, 0, 1,
                          rmatind.size(), &rhs, &sense, &rmatbeg,
                          &rmatind[0], &rmatval[0], (char **) NULL,
                          (char **) NULL);

    if (rval)
        throw cpx_err(rval, "CPXaddrows");
}

void Relaxation::add_cuts(const vector<double> &rhs,
                          const vector<char> &sense,
                          const vector<int> &rmatbeg,
                          const vector<int> &rmatind,
                          const vector<double> &rmatval)
{
    int rval = CPXaddrows(simpl_p->env, simpl_p->lp, 0, rmatbeg.size(),
                          rmatind.size(), &rhs[0], &sense[0],
                          &rmatbeg[0], &rmatind[0], &rmatval[0],
                          (char **) NULL, (char **) NULL);
    if (rval)
        throw cpx_err(rval, "CPXaddrows");
}

void Relaxation::del_set_rows(std::vector<int> &delstat)
{
    int rval = CPXdelsetrows(simpl_p->env, simpl_p->lp, &delstat[0]);
    if (rval)
        throw cpx_err(rval, "CPXdelsetrows");
}

void Relaxation::get_row_infeas(const std::vector<double> &x,
                                std::vector<double> &feas_stat,
                                int begin, int end) const
{
    feas_stat.resize(end - begin + 1);

    int rval = CPXgetrowinfeas(simpl_p->env, simpl_p->lp, &x[0], &feas_stat[0],
                               begin, end);

    if (rval)
        throw cpx_err(rval, "CPXgetrowinfeas");
}

void Relaxation::add_col(const double objval, const vector<int> &indices,
                         const vector<double> &coeffs, const double lb,
                         const double ub)
{
    int cmatbeg = 0;
    int newcols = 1;

    int rval = CPXaddcols(simpl_p->env, simpl_p->lp, newcols, coeffs.size(),
                          &objval, &cmatbeg, &indices[0], &coeffs[0],
                          &lb, &ub, (char **) NULL);
    if (rval)
        throw cpx_err(rval, "CPXaddcols");
}

void Relaxation::get_base(vector<int> &colstat,
                          vector<int> &rowstat) const
{
    colstat.resize(num_cols());
    rowstat.resize(num_rows());
    
    int rval = CPXgetbase(simpl_p->env, simpl_p->lp, &colstat[0], &rowstat[0]);
    if (rval)
        throw cpx_err(rval, "CPXgetbase");
}

Basis Relaxation::base() const
{
    Basis result;
    get_base(result.colstat, result.rowstat);
    return result;
}

vector<int> Relaxation::col_stat() const
{
    vector<int> result(num_cols());

    int rval = CPXgetbase(simpl_p->env, simpl_p->lp, &result[0], NULL);
    if (rval)
        throw cpx_err(rval, "CPXgetbase");

    return result;
}

void Relaxation::copy_start(const vector<double> &x)
{
    int rval = CPXcopystart(simpl_p->env, simpl_p->lp,
                            (int *) NULL, (int *) NULL, //starting basis
                            &x[0], (double *) NULL, //starting primal vals
                            (double *) NULL, (double *) NULL); //starting duals
    if (rval)
        throw cpx_err(rval, "CPXcopystart");
}

/**
 * @param[in] col_stat the basic statuses for the columns.
 * @param[in] row_stat the basic statuses for the rows.
 */
void Relaxation::copy_base(const std::vector<int> &col_stat,
                           const std::vector<int> &row_stat)
{
    int rval = CPXcopybase(simpl_p->env, simpl_p->lp, &col_stat[0],
                           &row_stat[0]);
    if (rval)
        throw cpx_err(rval, "CPXcopybase");
}

/**
 * @param[in] base the Basis structure containing the row and column statuses
 * to be copied.
 */
void Relaxation::copy_base(const Basis &base)
{
    copy_base(base.colstat, base.rowstat);
}

void Relaxation::copy_start(const vector<double> &x,
                            const vector<int> &col_stat,
                            const vector<int> &row_stat)
{
    int rval = CPXcopystart(simpl_p->env, simpl_p->lp,
                            &col_stat[0], &row_stat[0], //starting basis
                            &x[0], (double *) NULL, //starting primal vals
                            (double *) NULL, (double *) NULL); //starting duals
    if (rval)
        throw cpx_err(rval, "CPXcopystart");
}

void Relaxation::factor_basis()
{
    int rval = 0;
    CPXlongParamGuard it_clamp_guard(CPX_PARAM_ITLIM, 0, simpl_p->env,
                                     "factor basis itlim clamp");

    rval = CPXprimopt(simpl_p->env, simpl_p->lp);
    if (rval)
        throw cpx_err(rval, "CPXprimopt 0 iterations");

    int solstat = CPXgetstat(simpl_p->env, simpl_p->lp);
    if (solstat != CPX_STAT_ABORT_IT_LIM)
        throw cpx_err(solstat, "CPXgetstat in factor_basis");
}

void Relaxation::primal_opt()
{
    int rval = CPXprimopt(simpl_p->env, simpl_p->lp);
    if (rval)
        throw cpx_err(rval, "CPXprimopt");
}

void Relaxation::dual_opt()
{
    int rval = CPXdualopt(simpl_p->env, simpl_p->lp);
    if (rval)
        throw cpx_err(rval, "CPXdualopt");
}

void Relaxation::nondegen_pivot(const double lowlimit)
{
    runtime_error err("Problem in Relaxation::nondegen_pivot.");
    int rval = 0;

    CPXdblParamGuard obj_ll(CPX_PARAM_OBJLLIM, lowlimit, simpl_p->env,
                            "nondegen_pivot obj limit");

    rval = CPXprimopt(simpl_p->env, simpl_p->lp);
    if (rval) {
        cerr << "CPXprimopt failed, rval: " << rval << "\n";
        throw err;
    }

    int solstat = CPXgetstat(simpl_p->env, simpl_p->lp);
    if (solstat == CPX_STAT_INFEASIBLE) {
        cerr << "Relaxation is infeasible.\n";
        throw err;
    }

    if (solstat != CPX_STAT_OPTIMAL &&
        solstat != CPX_STAT_OPTIMAL_INFEAS &&
        solstat != CPX_STAT_ABORT_OBJ_LIM ) {
        cerr << "Solstat: " << solstat << "\n";
        throw err;
    }    
}

void Relaxation::primal_pivot()
{
    int rval = 0;
    CPXlongParamGuard it_clamp(CPX_PARAM_ITLIM, 1, simpl_p->env,
                               "single pivot itlim clamp");

    rval = CPXprimopt(simpl_p->env, simpl_p->lp);
    if (rval)
        throw cpx_err(rval, "CPXprimopt");

    int solstat = CPXgetstat(simpl_p->env, simpl_p->lp);
    if (solstat == CPX_STAT_INFEASIBLE)
        throw runtime_error("LP is infeasible.");

    if (solstat != CPX_STAT_OPTIMAL &&
        solstat != CPX_STAT_ABORT_IT_LIM &&
        solstat != CPX_STAT_OPTIMAL_INFEAS)
        throw cpx_err(solstat, "CPXprimopt solstat");
    
}

void Relaxation::dual_pivot()
{
    int rval = 0;
    CPXlongParamGuard it_clamp(CPX_PARAM_ITLIM, 1, simpl_p->env,
                               "single pivot itlim clamp");

    rval = CPXdualopt(simpl_p->env, simpl_p->lp);
    if (rval)
        throw cpx_err(rval, "CPXdualopt");

    int solstat = CPXgetstat(simpl_p->env, simpl_p->lp);
    if (solstat == CPX_STAT_INFEASIBLE)
        throw runtime_error("LP is infeasible.");

    if (solstat != CPX_STAT_OPTIMAL &&
        solstat != CPX_STAT_ABORT_IT_LIM &&
        solstat != CPX_STAT_OPTIMAL_INFEAS)
        throw cpx_err(solstat, "CPXdualopt solstat");
    
}

/**
 * Using the resident basis as a starting point, this function will invoke
 * the primal simplex optimizer, stopping just as soon as a pivot renders the
 * basis primal feasible.
 */
void Relaxation::primal_recover()
{
    int rval = CPXsetlpcallbackfunc(simpl_p->env, pfeas_cb, NULL);
    if (rval)
        throw cpx_err(rval, "CPXsetlpcallbackfunc setting pfeas_cb");

    primal_opt();

    rval = CPXsetlpcallbackfunc(simpl_p->env, NULL, NULL);
    if (rval)
        throw cpx_err(rval, "CPXsetlpcallbackfunc undoing cb");
}

double Relaxation::get_objval() const
{
    double result = std::numeric_limits<double>::max();
    
    int rval = CPXgetobjval(simpl_p->env, simpl_p->lp, &result);
    if (rval)
        throw cpx_err(rval, "CPXgetobjval");

    return result;    
}

void Relaxation::get_x(vector<double> &x) const
{

    set_info_vec(CPXgetx, "CPXgetx", simpl_p->env, simpl_p->lp, x, 0,
                 num_cols() - 1);
}

vector<double> Relaxation::lp_vec() const
{
    return info_vec(CPXgetx, "CPXgetx", simpl_p->env, simpl_p->lp, 0,
                    num_cols() - 1);
}

bool Relaxation::primal_feas() const
{
    int result = 0;
    int rval = CPXsolninfo(simpl_p->env, simpl_p->lp, NULL, NULL, &result,
                           NULL);
    if (rval)
        throw cpx_err(rval, "CPXsolninfo");
    return result;
}

bool Relaxation::dual_feas() const
{
    int result = 0;

    int rval = CPXsolninfo(simpl_p->env, simpl_p->lp, NULL, NULL, NULL,
                           &result);
    if (rval)
        throw cpx_err(rval, "CPXsolninfo");

    return result;
}

void Relaxation::get_row_slacks(vector<double> &slack, int begin,
                                int end) const
{
    set_info_vec(CPXgetslack, "CPXgetslack", simpl_p->env, simpl_p->lp,
                 slack, begin, end);
}

vector<double> Relaxation::row_slacks(int begin, int end) const
{
    return info_vec(CPXgetslack, "CPXgetslack", simpl_p->env, simpl_p->lp,
                    begin, end);
}

void Relaxation::get_pi(vector<double> &pi, int begin, int end) const
{
    set_info_vec(CPXgetpi, "CPXgetpi", simpl_p->env, simpl_p->lp, pi, begin,
                 end);
}

vector<double> Relaxation::pi(int begin, int end) const
{
    return info_vec(CPXgetpi, "CPXgetpi", simpl_p->env, simpl_p->lp, begin,
                    end);
}

void Relaxation::get_redcosts(vector<double> &rcs, int begin, int end) const
{
    set_info_vec(CPXgetdj, "CPXgetdj", simpl_p->env, simpl_p->lp, rcs, begin,
                 end);
}

vector<double> Relaxation::redcosts(int begin, int end) const
{
    return info_vec(CPXgetdj, "CPXgetdj", simpl_p->env, simpl_p->lp, begin,
                    end);
}

/**
 * @param[in] tour_vec the resident best tour
 * @param[in] colstat the column basis for \p tour_vec
 * @param[in] rowstat the row basis for \p tour_vec
 * @param[in] indices the column indices of edges to examine
 * @param[out] downobj estimates for clamping a variable to zero.
 * @param[out] upobj estimates for clamping a variable to one. 
 * @param[in/out] contra_bases primal feasible starting bases for enforcing 
 * Contra branches. May be empty, in which case it will be populated for use
 * in future calls. If nonempty, it will be assumed to have the same size as
 * \p indices, with `contra_bases[i]` to be used as a starting basis for fixing
 * `indices[i]` to disagree with `tour_entry[indices[i]]`.
 * @param[in] itlim the maximum number of simplex iterations to do.
 * @param[in] upperbound the length of \p tour_vec.
 */
void Relaxation::primal_strong_branch(const vector<double> &tour_vec,
                                      const vector<int> &colstat,
                                      const vector<int> &rowstat,
                                      const vector<int> &indices,
                                      vector<std::pair<int, double>> &downobj,
                                      vector<std::pair<int, double>> &upobj,
                                      vector<Basis> &contra_bases,
                                      int itlim, double upperbound)
{
    using ScorePair = std::pair<int, double>;
    
    CPXintParamGuard per_ind(CPX_PARAM_PERIND, 0, simpl_p->env,
                             "primal_strong_branch perturb");

    CPXintParamGuard price_ind(CPX_PARAM_PPRIIND, CPX_PPRIIND_STEEP,
                               simpl_p->env, "primal_strong_branch pricing");
    
    downobj.clear();
    upobj.clear();
    downobj.reserve(indices.size());
    upobj.reserve(indices.size());

    bool have_bases = false;
    if (!contra_bases.empty())
        have_bases = true;
    else
        contra_bases.reserve(indices.size());

    using ClampPair = std::pair<char, double>;

    std::array<ClampPair, 2> clamps{ClampPair('U', 0.0), ClampPair('L', 1.0)};

    for (int i = 0; i < indices.size(); ++i) {
        int ind = indices[i];
        //cout << "----Index " << ind << "\n";
        for (ClampPair &cp : clamps) {
            char sense = cp.first;
            double clamp_bound = cp.second;
            double unclamp_bound = 1.0 - cp.second;

            tighten_bound(ind, sense, clamp_bound);

            if (tour_vec[ind] == clamp_bound) {
                copy_start(tour_vec);
                factor_basis();
            } else {
                if (!have_bases) {
                    copy_base(colstat, rowstat);
                    primal_recover();
                    // cout << "P feas after prim recover: "
                    //      << primal_feas() << ", "
                    //      << it_count() << " iterations\n";
                    if (!primal_feas())
                        cout << "Infeasible with stat "
                             << CPXgetstat(simpl_p->env, simpl_p->lp) << "\n";
                    contra_bases.emplace_back(base());
                } else {
                    copy_base(contra_bases[i].colstat,
                              contra_bases[i].rowstat);
                    factor_basis();
                    // cout << "Copied saved base, p feas "
                    //      << primal_feas() << "\n";
                }
            }

            CPXlongParamGuard it_lim(CPX_PARAM_ITLIM, itlim, simpl_p->env,
                                     "primal_strong_branch it lim");

            // cout << "objval after tighten " << ((int) clamp_bound) << ": "
            //      << get_objval() << "\n\n";

            primal_opt();

            int solstat = CPXgetstat(simpl_p->env, simpl_p->lp);
            double objval = get_objval();
            int rank = -1;

            if (solstat == CPX_STAT_ABORT_IT_LIM)
                rank = 0;
            else if (solstat == CPX_STAT_OPTIMAL ||
                     solstat == CPX_STAT_OPTIMAL_INFEAS)
                rank = 1;
            else if (solstat == CPX_STAT_INFEASIBLE)
                rank = 2;
            else 
                throw cpx_err(solstat,
                              clamp_bound == 0.0 ?
                              "CPXgetstat in down clamp" :
                              "CPXgetstat in up clamp");

            vector<ScorePair> &objvec = clamp_bound == 0.0 ? downobj : upobj;
            objvec.push_back(ScorePair(rank, objval));

            tighten_bound(ind, sense, unclamp_bound);
        }
    }

    copy_start(tour_vec);
    factor_basis();
}

void Relaxation::tighten_bound(const int index, const char sense,
                               const double val)
{
    int rval = CPXtightenbds(simpl_p->env, simpl_p->lp, 1, &index, &sense,
                             &val);
    if (rval)
        throw cpx_err(rval, "CPXtightenbds");
}

void Relaxation::change_obj(const int index, const double val)
{
    int rval = CPXchgobj(simpl_p->env, simpl_p->lp, 1, &index, &val);
    if (rval)
        throw cpx_err(rval, "CPXchgobj");
}

#if CMR_HAVE_SAFEGMI

void Relaxation::init_mir_data(Sep::MIRgroup &mir_data)
{
    runtime_error err("Problem in Relaxation::init_mir_data.");

    using mir_lp = SLVRcplex_t;
    using mir_system = CUTSsystem_t<double>;
    using mir_varinfo = CUTSvarInfo_t<double>;
    using mir_basinfo = SLVRbasisInfo_t;

    // Initialize the solver/constraint info //
    
    int numcols = num_cols();
    vector<char> ctype;

    try {
        ctype = vector<char>(numcols, 'B');
    } CMR_CATCH_PRINT_THROW("allocating ctype", err);

    mir_lp lp_obj;
    lp_obj.env = simpl_p->env;
    lp_obj.prob = simpl_p->lp;
    lp_obj.ctype = &ctype[0];

    mir_system *constraint_matrix = (mir_system *) NULL;

    if(SLVRformulationRows(&lp_obj, &constraint_matrix)) {
        cerr << "SLVRformulationrows failed.\n";
        throw err;
    }

    mir_data.constraint_matrix =
    unique_ptr<mir_system, Sep::SystemDeleter>(constraint_matrix);

    mir_basinfo *binfo = (mir_basinfo *) NULL;
    if (SLVRgetBasisInfo(&lp_obj, &binfo)) {
        cerr << "SLVRgetBasisInfo failed.\n";
        throw err;
    }

    struct BinfoDeleter {
        void operator()(mir_basinfo *P) {
            SLVRfreeBasisInfo(&P);
        }
    };

    unique_ptr<mir_basinfo, BinfoDeleter> basis_info(binfo);

    mir_varinfo *vinfo = (mir_varinfo *) NULL;
    if (SLVRgetVarInfo(&lp_obj, true, &vinfo)) {
        cerr << "SLVRgetVarInfo failed.\n";
        throw err;
    }

    mir_data.var_info = unique_ptr<mir_varinfo, Sep::VinfoDeleter>(vinfo);

    vector<double> x;
    try { get_x(x); } CMR_CATCH_PRINT_THROW("getting x", err);

    mir_data.full_x =
    util::c_array_ptr<double>(SLVRgetFullX(&lp_obj,
                                           mir_data.constraint_matrix.get(),
                                           &x[0]));
    
    if (mir_data.full_x.get() == NULL) {
        cerr << "SLVRgetFullX failed.\n";
        throw err;
    }

    // now rank fractional basic structural variables and get tableau rows //

    using VarPair = std::pair<int, double>;
    vector<double> &lp_vranking = mir_data.vranking;
    vector<int> colstat;
    vector<VarPair> frac_basic_vars;

    try {
        //only want cuts from structural variables but the vector needs to be
        //size of full_x, so those are just left with a -1 ranking.
        lp_vranking.resize(num_rows() + num_cols(), -1.0);
        colstat = col_stat();
    } CMR_CATCH_PRINT_THROW("resizing and getting col stat", err);

    for (int i = 0; i < x.size(); ++i) {
        if (colstat[i] == 1 && !util::var_integral(x[i])) {
            lp_vranking[i] = (-(x[i] - 0.5) * (x[i] - 0.5)) + 0.25;
            try {
                frac_basic_vars.emplace_back(i, x[i]);
            } CMR_CATCH_PRINT_THROW("pushing back frac basic var", err);
        }
    }

    if (frac_basic_vars.empty())
        throw logic_error("Tried init_mir_data w no fractional basic vars.");

    std::sort(frac_basic_vars.begin(), frac_basic_vars.end(),
              [](VarPair a, VarPair b)
              { return abs(0.5 - a.second) < abs(0.5 - b.second); });

    mir_system *tab_rows = (mir_system *) NULL;
    if (CUTSnewSystem(&tab_rows, frac_basic_vars.size())) {
        cerr << "CUTSnewSystem failed.\n";
        throw err;
    }

    for (VarPair v : frac_basic_vars) {
        if (SLVRgetTableauRow(&lp_obj,
                              &constraint_matrix,
                              &(tab_rows->rows[tab_rows->sys_rows]),
                              &binfo,
                              v.first)) {
            CUTSfreeSystem(&tab_rows);
            cerr << "SLVRgetTableauRow failed.\n";
            throw err;
        }

        tab_rows->sys_rows += 1;
    }

    mir_data.tableau_rows = unique_ptr<mir_system,
                                       Sep::SystemDeleter>(tab_rows);
}

#endif //CMR_HAVE_SAFEGMI


}
}
