#include "lp_interface.hpp"
#include "err_util.hpp"

#include <stdexcept>
#include <iostream>
#include <limits>

using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::vector;

using cpx_err = CMR::retcode_error;


namespace CMR {
namespace LP {

Relaxation::Relaxation() try
{
    int rval = 0;
    
    cplex_lp = (CPXLPptr) NULL;
    cplex_env = CPXopenCPLEX(&rval);

    if (rval) 
        throw cpx_err(rval, "CPXopenCPLEX");

    
    rval = CPXsetintparam(cplex_env, CPX_PARAM_PREIND, CPX_OFF);

    if (rval)
        throw cpx_err(rval, "CPXsetintparam presolve");

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

int Relaxation::num_rows() { return CPXgetnumrows(cplex_env, cplex_lp); }

int Relaxation::num_cols() { return CPXgetnumcols(cplex_env, cplex_lp); }

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

void Relaxation::get_row_infeas(const std::vector<double> &x,
                                std::vector<double> &feas_stat,
                                int begin, int end)
{
    feas_stat.resize(end - begin + 1);

    int rval = CPXgetrowinfeas(cplex_env, cplex_lp, &x[0], &feas_stat[0],
                               begin, end);

    if (rval)
        throw cpx_err(rval, "CPXgetrowinfeas");
}

void Relaxation::new_cols(const vector<double> &objfun,
                          const vector<int> &cmatbeg,
                          const vector<int> &cmatind,
                          const vector<double> &cmatval,
                          const vector<double> &lb, const vector<double> &ub)
{
    int rval = CPXnewcols(cplex_env, cplex_lp, cmatbeg.size(),
                          &objfun[0], &lb[0], &ub[0], (char *) NULL,
                          (char **) NULL);
    if (rval)
        throw cpx_err(rval, "CPXnewcols");
}

void Relaxation::get_base(vector<int> &colstat, vector<int> &rowstat)
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

double Relaxation::get_objval()
{
    double result = std::numeric_limits<double>::max();
    
    int rval = CPXgetobjval(cplex_env, cplex_lp, &result);
    if (rval)
        throw cpx_err(rval, "CPXgetobjval");

    return result;    
}

void Relaxation::get_x(vector<double> &x)
{
    x.resize(num_cols());
    
    int rval = CPXgetx(cplex_env, cplex_lp, &x[0], 0, num_cols() - 1);
    if (rval)
        throw cpx_err(rval, "CPXgetx");
}

bool Relaxation::dual_feas()
{
    int result = 0;

    int rval = CPXsolninfo(cplex_env, cplex_lp, NULL, NULL, NULL,
                           &result);
    if (rval)
        throw cpx_err(rval, "CPXsolninfo");

    return result;
}

void Relaxation::get_row_slacks(vector<double> &slack, int begin, int end)
{
    slack.resize(end - begin + 1);

    int rval = CPXgetslack(cplex_env, cplex_lp, &slack[0], begin, end);
    if (rval)
        throw cpx_err(rval, "CPXgetslack");
}

void Relaxation::get_pi(vector<double> &pi, int begin, int end)
{
    pi.resize(end - begin + 1);

    int rval = CPXgetpi(cplex_env, cplex_lp, &pi[0], begin, end);
    if (rval)
        throw cpx_err(rval, "CPXgetpi");
}

}
}
