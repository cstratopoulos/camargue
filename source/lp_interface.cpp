#include "lp_interface.hpp"
#include "err_util.hpp"

#include <stdexcept>
#include <iostream>
#include <limits>
#include <typeinfo>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::to_string;

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::vector;

using cpx_err = CMR::retcode_error;


namespace CMR {
namespace LP {

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

template <typename cplex_query>
void set_info_vec(cplex_query F, const char *Fname,
                  const CPXENVptr cplex_env, const CPXLPptr cplex_lp,
                  vector<double> &info_vec, int begin, int end)
{
    info_vec.resize(end - begin + 1);

    int rval = F(cplex_env, cplex_lp, &info_vec[0], begin, end);
    if (rval)
        throw cpx_err(rval, Fname);
}

Relaxation::Relaxation() try
{
    int rval = 0;
    
    cplex_lp = (CPXLPptr) NULL;
    cplex_env = CPXopenCPLEX(&rval);

    if (rval) 
        throw cpx_err(rval, "CPXopenCPLEX");

    rval = CPXsetintparam(cplex_env, CPX_PARAM_AGGIND, CPX_OFF);
    if (rval) {
        CPXcloseCPLEX(&cplex_env);
        throw cpx_err(rval, "CPXsetintparam aggregator");
    }
    
    rval = CPXsetintparam(cplex_env, CPX_PARAM_PREIND, CPX_OFF);
    if (rval) {
        CPXcloseCPLEX(&cplex_env);
        throw cpx_err(rval, "CPXsetintparam presolve");
    }

    rval = CPXsetintparam(cplex_env, CPX_PARAM_PPRIIND, CPX_PPRIIND_DEVEX);
    if (rval) {
        CPXcloseCPLEX(&cplex_env);
        throw cpx_err(rval, "CPXsetintparam primal pricing");
    }

    
    char unused;

    cplex_lp = CPXcreateprob(cplex_env, &rval, &unused);

    if (rval) {
        CPXcloseCPLEX(&cplex_env);
        throw cpx_err(rval, "CPXcreateprob");
    }

    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Relaxation constructor failed.");
}

Relaxation::Relaxation(Relaxation &&lp) noexcept : cplex_env(lp.cplex_env),
                                       cplex_lp(lp.cplex_lp)
{
    if (cplex_env) {
        if (cplex_lp) {
            CPXfreeprob(cplex_env, &cplex_lp);
            cplex_lp = (CPXLPptr) NULL;
        }
        CPXcloseCPLEX(&cplex_env);
        cplex_env = (CPXENVptr) NULL;
    }    
    
    cplex_env = lp.cplex_env;
    cplex_lp = lp.cplex_lp;

    lp.cplex_env = (CPXENVptr) NULL;
    lp.cplex_lp = (CPXLPptr) NULL;
}

Relaxation& Relaxation::operator=(Relaxation &&lp) noexcept
{
    if (cplex_env) {
        if (cplex_lp) {
            CPXfreeprob(cplex_env, &cplex_lp);
            cplex_lp = (CPXLPptr) NULL;
        }
        CPXcloseCPLEX(&cplex_env);
        cplex_env = (CPXENVptr) NULL;
    }
    
    cplex_env = lp.cplex_env;
    cplex_lp = lp.cplex_lp;
    lp.cplex_env = (CPXENVptr) NULL;
    lp.cplex_lp = (CPXLPptr) NULL;
    return *this;
}

Relaxation::~Relaxation()
{
    if (cplex_env) {
        if (cplex_lp) {
            CPXfreeprob(cplex_env, &cplex_lp);
            cplex_lp = (CPXLPptr) NULL;
        }
        CPXcloseCPLEX(&cplex_env);
        cplex_env = (CPXENVptr) NULL;
    }
}

int Relaxation::num_rows() const { return CPXgetnumrows(cplex_env, cplex_lp); }

int Relaxation::num_cols() const { return CPXgetnumcols(cplex_env, cplex_lp); }

void Relaxation::new_row(const char sense, const double rhs)
{
    int rval = CPXnewrows(cplex_env, cplex_lp, 1, &rhs, &sense, NULL, NULL);

    if (rval)
        throw cpx_err(rval, "CPXnewrows");
}

void Relaxation::new_rows(const vector<char> &sense,
                          const vector<double> &rhs)
{
    int rval = CPXnewrows(cplex_env, cplex_lp, sense.size(), rhs.data(),
                          sense.data(), NULL, NULL);
    if (rval)
        throw cpx_err(rval, "CPXnewrows");
}

void Relaxation::add_cut(const double rhs, const char sense,
                         const vector<int> &rmatind,
                         const vector<double> &rmatval)
{
    int rmatbeg = 0;

    int rval = CPXaddrows(cplex_env, cplex_lp, 0, 1,
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
    int rval = CPXaddrows(cplex_env, cplex_lp, 0, rmatbeg.size(),
                          rmatind.size(), &rhs[0], &sense[0],
                          &rmatbeg[0], &rmatind[0], &rmatval[0],
                          (char **) NULL, (char **) NULL);
    if (rval)
        throw cpx_err(rval, "CPXaddrows");
}

void Relaxation::del_set_rows(std::vector<int> &delstat)
{
    int rval = CPXdelsetrows(cplex_env, cplex_lp, &delstat[0]);
    if (rval)
        throw cpx_err(rval, "CPXdelsetrows");
}

void Relaxation::get_row_infeas(const std::vector<double> &x,
                                std::vector<double> &feas_stat,
                                int begin, int end) const
{
    feas_stat.resize(end - begin + 1);

    int rval = CPXgetrowinfeas(cplex_env, cplex_lp, &x[0], &feas_stat[0],
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

    int rval = CPXaddcols(cplex_env, cplex_lp, newcols, coeffs.size(),
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
    
    int rval = CPXgetbase(cplex_env, cplex_lp, &colstat[0], &rowstat[0]);
    if (rval)
        throw cpx_err(rval, "CPXgetbase");
}

void Relaxation::copy_start(const vector<double> &x)
{
    int rval = CPXcopystart(cplex_env, cplex_lp, (int *) NULL, (int *) NULL,
                            &x[0], (double *) NULL, (double *) NULL,
                            (double *) NULL);
    if (rval)
        throw cpx_err(rval, "CPXcopystart");
}

void Relaxation::copy_start(const vector<double> &x,
                            const vector<int> &col_stat,
                            const vector<int> &row_stat)
{
    int rval = CPXcopystart(cplex_env, cplex_lp, &col_stat[0],
                            &row_stat[0], &x[0], (double *) NULL,
                            (double *) NULL, (double *) NULL);
    if (rval)
        throw cpx_err(rval, "CPXcopystart");
}

void Relaxation::factor_basis()
{
    runtime_error err("problem in factor_basis.");

    int rval = CPXsetlongparam(cplex_env, CPX_PARAM_ITLIM, 0);
    if (rval) {
        cerr << "CPXsetlongparam itlim clamp failed, rval: " << rval << "\n";
        throw err;
    }

    rval = CPXprimopt(cplex_env, cplex_lp);
    if (rval) {
        cerr << "CPXprimopt failed, rval: " << rval << "\n";
        throw err;
    }

    int solstat = CPXgetstat(cplex_env, cplex_lp);
    if (solstat != CPX_STAT_ABORT_IT_LIM) {
        cerr << "Solstat: " << solstat << "\n";
        throw err;
    }

    rval = CPXsetlongparam(cplex_env, CPX_PARAM_ITLIM,
                           9223372036800000000);
    if (rval) {
        cerr << "CPXsetlong param itlim revert failed, rval: " << rval << "\n";
        throw err;
    }
}

void Relaxation::primal_opt()
{
    int rval = CPXprimopt(cplex_env, cplex_lp);
    if (rval)
        throw cpx_err(rval, "CPXprimopt");
}

void Relaxation::nondegen_pivot(const double lowlimit)
{
    runtime_error err("Problem in Relaxation::nondegen_pivot.");
    
    int rval = CPXsetdblparam(cplex_env, CPX_PARAM_OBJLLIM, lowlimit);
    
    if (rval) {
        cerr << "CPXsetdblparam failed setting low limit, rval: "
             << rval << "\n";
        throw err;
    }

    rval = CPXprimopt(cplex_env, cplex_lp);
    if (rval) {
        cerr << "CPXprimopt failed, rval: " << rval << "\n";
        throw err;
    }

    int solstat = CPXgetstat(cplex_env, cplex_lp);
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

    rval = CPXsetdblparam(cplex_env, CPX_PARAM_OBJLLIM, -1E75);
    if (rval) {
        cerr << "CPXsetdblparam failed reverting low limit, rval: "
             << rval << "\n";
        throw err;
    }    
}

double Relaxation::get_objval() const
{
    double result = std::numeric_limits<double>::max();
    
    int rval = CPXgetobjval(cplex_env, cplex_lp, &result);
    if (rval)
        throw cpx_err(rval, "CPXgetobjval");

    return result;    
}

void Relaxation::get_x(vector<double> &x) const
{

    set_info_vec(CPXgetx, "CPXgetx", cplex_env, cplex_lp, x, 0,
                 num_cols() - 1);
}

vector<double> Relaxation::lp_vec() const
{
    return info_vec(CPXgetx, "CPXgetx", cplex_env, cplex_lp, 0,
                    num_cols() - 1);
}


bool Relaxation::dual_feas() const
{
    int result = 0;

    int rval = CPXsolninfo(cplex_env, cplex_lp, NULL, NULL, NULL,
                           &result);
    if (rval)
        throw cpx_err(rval, "CPXsolninfo");

    return result;
}

void Relaxation::get_row_slacks(vector<double> &slack, int begin,
                                int end) const
{
    set_info_vec(CPXgetslack, "CPXgetslack", cplex_env, cplex_lp,
                 slack, begin, end);
}

vector<double> Relaxation::row_slacks(int begin, int end) const
{
    return info_vec(CPXgetslack, "CPXgetslack", cplex_env, cplex_lp,
                    begin, end);
}

void Relaxation::get_pi(vector<double> &pi, int begin, int end) const
{
    set_info_vec(CPXgetpi, "CPXgetpi", cplex_env, cplex_lp, pi, begin, end);
}

vector<double> Relaxation::pi(int begin, int end) const
{
    return info_vec(CPXgetpi, "CPXgetpi", cplex_env, cplex_lp, begin, end);
}

}
}
