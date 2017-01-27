/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Configuration macros.
 *
 * This header file contains macro definitions that are used to conditionally
 * compile sections of the code based on presence or absence of certain 
 * external dependences.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CONFIG_H
#define CMR_CONFIG_H

/** @name Existence macros
 * These govern indicate certain features are downloaded or supported by the 
 * compiler.
 */
///@{

/// Define if you have downloaded the header file for Catch unit testing.
/// Catch is the unit testing framework used in this project. If you would like
/// to compile and run this project's unit tests yourself, Catch can be
/// obtained at https://github.com/philsquared/Catch
#define CMR_HAVE_CATCH 1


/// Define if you have downloaded the header for Timsort implementation.
/// As a heuristic enhancement, this project uses the Timsort sorting algorithm
/// to perform certain sorting operations in generating candidate teeth. If
/// you would like to use Timsort in these cases, an implementation by
/// Fuji Goro can be obtained at https://github.com/gfx/cpp-TimSort
#define CMR_HAVE_TIMSORT 1

/// Define if your compiler supports OpenMP.
/// Some sections of this project implement simple algorithms/subroutines in
/// parallel by use of OpenMP (OMP). The OMP specification guarantees that if
/// the compiler does not support OMP, the pragmas are simply ignored by the
/// compiler and the result is still valid code. However, in some cases there
/// is a bit of added overhead for error handling or memory allocation which
/// is not necessary in the serial case.
#define CMR_HAVE_OPENMP 1

#define CMR_HAVE_SAFEGMI 1

///@}



/** @name Usage macros
 * Timsort is always used if it is present, but more care is required with OMP
 * and Catch: Catch is not threadsafe, so running Catch tests with OMP enabled
 * can cause memory leaks with thread creation/destruction. Thus, OMP is used
 * iff the compiler supports it and tests are not being ran.
 */
///@{

#ifdef CMR_HAVE_CATCH

/// Define if you want to compile tests and invoke them from the command line.
//#define CMR_DO_TESTS

#endif //CMR_HAVE_CATCH



#ifndef CMR_DO_TESTS
#ifdef CMR_HAVE_OPENMP

/// Define if OMP parallelism will actually be used in the code. 
#define CMR_USE_OMP

#endif //CMR_HAVE_OPENMP
#endif //CMR_DO_TESTS


///@}




#endif // CMR_CONFIG_H
