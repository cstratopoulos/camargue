/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Utilities for working with CPLEX.
 * Class and function templates in this header generally provide more C++-like
 * semantics and interfaces for the CPLEX Callable library.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CPX_UTIL_H
#define CMR_CPX_UTIL_H

#include "err_util.hpp"

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>

#include <cplex.h>

namespace CMR {
namespace LP {

using cpx_err = CMR::util::retcode_error;

/** Template for getting a ranged result vector from CPLEX.
 * @tparam cplex_query the function type of the query.
 * @param F the function to call.
 * @param Fname the function name, if an error message needs to be printed.
 * @returns a vector obtained by calling \p F to obtain info for the range
 * \p begin to \p end.
 */
template <typename cplex_query>
std::vector<double> info_vec(cplex_query F, const char* Fname,
                             CPXENVptr cplex_env, CPXLPptr cplex_lp,
                             int begin, int end)
{
    std::vector<double> result(end - begin + 1);

    int rval = F(cplex_env, cplex_lp, &result[0], begin, end);
    if (rval)
        throw cpx_err(rval, Fname);

    return result;
}

/// Like info_vec, but modifies a vector rather than returning one.
template <typename cplex_query, typename vectype>
void set_info_vec(cplex_query F, const char *Fname,
                  CPXENVptr cplex_env, CPXLPptr cplex_lp,
                  std::vector<vectype> &info_vec, int begin, int end)
{
    info_vec.resize(end - begin + 1);

    int rval = F(cplex_env, cplex_lp, &info_vec[0], begin, end);
    if (rval)
        throw cpx_err(rval, Fname);
}

/** Gets callback info of a particular type from a user-written callback.
 * @tparam info_type the data type of the info to be returned.
 * See CPXgetcallbackinfo documentation for more info.
 */
template <typename info_type>
void get_callback_info(CPXCENVptr cpx_env, void *cbdata, int wherefrom,
                       int which_info, info_type *result_p,
                       const char *query_description)
{
    int rval = CPXgetcallbackinfo(cpx_env, cbdata, wherefrom, which_info,
                                  result_p);
    if (rval)
        throw cpx_err(rval, query_description);
}

/// Optimization callback for Relaxation::primal_recover.
static int pfeas_cb(CPXCENVptr cpx_env, void *cbdata, int wherefrom,
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

/// implementation aliases/details for cplex param guards.
namespace detail {

/// Template alias for signature of CPXgetintparam, CPXgetdblparam, etc.
template<typename numtype>
using CPXgetType = int(*)(CPXCENVptr, int, numtype *);

/// Template alias for signature of CPXsetintparam, CPXsetdblparam, etc.
template<typename numtype>
using CPXsetType = int(*)(CPXENVptr, int, numtype);

/** A scope guard for making temporary changes to a CPLEX parameter.
 * @tparam numtype the numeric type of the parameter to change.
 * @tparam GetP a pointer to a function for getting a parameter of type
 * \p numtype.
 * @tparam SetP a pointer to a function for setting a parameter of type
 * \p numtype.
 */
template
<typename numtype, CPXgetType<numtype> GetP, CPXsetType<numtype> SetP>
class CPXparamGuard {
public:
    /// Construct a parameter guard. See private members for arguments.
    CPXparamGuard(int which, numtype new_value, CPXENVptr env,
                  const std::string p_desc) try
        : which_param(which), cplex_env(env), param_desc(p_desc)
    {
        int rval = GetP(cplex_env, which_param, &old_value);
        if (rval)
            throw cpx_err(rval, "Get param " + param_desc);

        rval = SetP(cplex_env, which_param, new_value);
        if (rval)
            throw cpx_err(rval, "Set param " + param_desc);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("CPXparamGuard constructor failed");
    }

    /// Revert the parameter, terminating the program upon failure.
    ~CPXparamGuard()
    {
        int rval = SetP(cplex_env, which_param, old_value);
        if (rval) {
            std::cerr << "\tFATAL: Failed to revert " << param_desc << ", rval "
                      << rval << " in destructor" << std::endl;
            exit(1);
        }
    }

private:
    const int which_param; //!< The CPLEX enumeration index of the parameter.
    CPXENVptr cplex_env; //!< The environemnt to change the param in.
    const std::string param_desc; //!< Description of the change being made.
    numtype old_value; //!< Set by constructor to the old value to revert to.
};
}

/// Integer parameter guard.
using CPXintParamGuard = detail::CPXparamGuard<int,
                                               &CPXgetintparam,
                                               &CPXsetintparam>;
/// Double parameter guard.
using CPXdblParamGuard = detail::CPXparamGuard<double,
                                               &CPXgetdblparam,
                                               &CPXsetdblparam>;
/// CPXLONG parameter guard.
using CPXlongParamGuard = detail::CPXparamGuard<CPXLONG,
                                                &CPXgetlongparam,
                                                &CPXsetlongparam>;


}
}

#endif
