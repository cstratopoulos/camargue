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

/** Macro for handling errors in function with multiple failure points.
 * Suggested usage: This macro should be used in a function by initially 
 * declaring a runtime_error err("Function failed."). Then, for various smaller
 * tasks within the function that may throw for distinct reasons (and which
 * throw exceptions derived from std::exception), use a try block followed by
 * CMR_CATCH_PRINT_THROW("description of failure in small task", err).
 */
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

    /// Construct a ScopeGuard to performs \p action upon going out of scope.
    ScopeGuard(act_type action) :
        final_action(std::move(action)), active(true) {}

    
    ~ScopeGuard()
        { if(active) final_action(); } //!< Perform the final_action.

    void dismiss()
        { active = false; } //!< Indicate that final_action should not happen.

    ScopeGuard() = delete; //!< No default constructor.
    ScopeGuard(const ScopeGuard&) = delete; //!< No copy constructor.
    ScopeGuard& operator=(const ScopeGuard&) = delete; //!< No copy assign.

    /// Move constructor. Moved from ScopeGuard \p rhs is dismissed.
    ScopeGuard(ScopeGuard&& rhs) :
        final_action(std::move(rhs.final_action)),
        active(rhs.active) {rhs.dismiss(); }
  
private:
    act_type final_action; //!< A function to be called on destruction.
    bool active; //!< True iff final_action should be performed on destruction.
};

/** Type deduction function for ScopeGuard, also from Andrei Alexandrescu.
 * @param[in] f the function to be called when the guard goes out of scope.
 * The suggested usage is `auto cleanup = make_guard([] {cleanup_lambda}); `
 */
template <typename act_type>
ScopeGuard<act_type> make_guard(act_type f)
{
  return ScopeGuard<act_type>(std::move(f));
}


/** Structure for converting retcodes to exceptions.
 * To be used with retcode-returning functions where different values signify
 * different errors, meaningful for lookup.
 */
struct retcode_error : public std::runtime_error {
    /**Construct a retcode_error with retcode and function name.
     * @param[in] rval the retcode from the function
     * @param[in] func_name the name of the function that may fail, or a brief
     * description of the task it performs.
     */
    retcode_error(const int rval, const std::string &func_name) :
        std::runtime_error(func_name + " failed with rval " +
                           std::to_string(rval)) {}
};

}  
}


#endif
