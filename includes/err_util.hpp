/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief ERROR HANDLING CLASSES AND ROUTINES
 *
 * This file contains classes for creating cleanup objects and exceptions 
 * derived from standard library runtime error. 
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_ERRUTIL_HPP
#define CMR_ERRUTIL_HPP

#include <stdexcept>
#include <string>
#include <iostream>

#define CMR_CATCH_PRINT_THROW( msg, new_ex )  catch(const std::exception &e)  \
  { std::cerr << e.what() << " " << msg << ".\n"; throw new_ex; }


namespace CMR {
namespace util {

/** Code from Andrei Alexandrescu's ScopeGuard11 slides.
 * This class provides a limited implementation of ScopeGuard. I have not 
 * bothered with the macro definitions for scope exit, or for operator + on 
 * a scope guard; all cleanup tasks must be specified upon construction. 
 */
template <typename act_type>
class ScopeGuard {
public:
  ScopeGuard(act_type action) :
    final_action(std::move(action)), active(true) {}
  ~ScopeGuard() { if(active) final_action(); }

  void dismiss() { active = false; }

  ScopeGuard() = delete;
  ScopeGuard(const ScopeGuard&) = delete;
  ScopeGuard& operator=(const ScopeGuard&) = delete;
  ScopeGuard(ScopeGuard&& rhs) :
    final_action(std::move(rhs.final_action)),
    active(rhs.active) {rhs.dismiss(); }
  
private:
  act_type final_action;
  bool active;
};

/** Type deduction function for ScopeGuard, also from Andrei Alexandrescu.
 * The suggested usage is `auto cleanup = make_guard([] {cleanup_lambda}); `
 */
template <typename act_type>
ScopeGuard<act_type> make_guard(act_type f)
{
  return ScopeGuard<act_type>(std::move(f));
}


/// Structure for converting retcodes to exceptions.
/// To be used with retcode-returning functions where different values signify
/// different errors, meaningful for lookup.
struct retcode_error : public std::runtime_error {
    /// Construct a retcode_error with retcode and function name.
    retcode_error(const int rval, const std::string &func_name) :
        std::runtime_error(func_name + " failed with rval " +
                           std::to_string(rval)) {}
};

}  
}


#endif
