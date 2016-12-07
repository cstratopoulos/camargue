/** @file
 *  @brief Header file for use of Catch tests.
 * This file contains a preprocessor definition used for conditional 
 * compilation of Catch tests. All files implementing test cases shall include
 * this file, and the catch.hpp header. 
 */


#ifndef CMR_TESTS_HPP
#define CMR_TESTS_HPP

/** Preprocessor definition for compiling with test cases.
 * If defined, compiling the project will create a Catch main for the
 * execution of test cases, rather than the main for running the solver. 
 */
#define CMR_DO_TESTS

#ifndef CMR_DO_TESTS

/** Preprocessor definition for use of OpenMP parallelism.
 * The use of this macro is twofold. 
 * 1. The OpenMP (OMP) specification guarantees that if the compiler does not
 * support OMP, then the pragmas are simply ignored and the result is still
 * valid code. However with some functions there is some extra overhead
 * related to  error handling or memory allocation that only occurs in the
 * parallel version, and it is possible to write a cleaner serial version, so
 * this one should be preferred.
 * 2. Camargue uses the Catch unit testing framework, which is not currently
 * threadsafe. On some tests, valgrind will show leaks related to thread
 * creation/deletion if parallelism is used. Thus, this macro's definition is
 * contingent on CMR_DO_TESTS being undefined. 
 */
#define CMR_USE_OMP

#endif


#endif
