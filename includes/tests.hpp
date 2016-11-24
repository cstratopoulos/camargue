/** @file
 *  @brief Header file for use of Catch tests.
 * This file contains a preprocessor definition used for conditional 
 * compilation of Catch tests. All files implementing test cases shall include
 * this file, and the catch.hpp header. 
 */


#ifndef PSEP_TESTS_HPP
#define PSEP_TESTS_HPP

/** Preprocessor definition for compiling with test cases.
 * If defined, compiling the project will create a Catch main for the
 * execution of test cases, rather than the main for running the solver. 
 */
//#define PSEP_DO_TESTS


#endif
