/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *               UTILITY FUNCTIONS AND MACROS/ENUM CLASSES
 *
 * This header contains namespace and enum classes for labels/constants that
 * are used elsewhere in the code, as well as some very simple structures
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

#include<utility>
#include<unordered_map>
#include<vector>
#include<iostream>
#include<type_traits>

#include<boost/functional/hash.hpp>

/*
 * In this project virtually all functions will produce a pseudo function
 * call graph in the event of an error, and they use the goto label CLEANUP
 * as the exit point for this. Also almost all functions have integer return 
 * types with nonzero values indicating error or special conditions. Thus
 * if the function big_task() contains a call to function helper_task(), 
 * use of this macro might look like
 *
 * int big_task()
 *    int rval = helper_task()
 *    if(rval) PSEP_GOTO_CLEANUP("helper_task failed, ");
 *    ... (skipped control flow) ...
 *    CLEANUP:
 *        if(rval) cerr << " problem in big_task\n"
 *
 * This will result in "helper_task failed, problem in big_task" being printed
 * to stderr, as well as similar error messages possibly embedded in 
 * helper_task, hopefully giving a clean description of the source of error
 */
#define PSEP_GOTO_CLEANUP(message) {std::cerr << message; goto CLEANUP;}
#define PSEP_CHECK_RVAL(rval, message) {\
  if ((rval)) {\
  std::cerr << message;\
  goto CLEANUP;\
  }\
  }

namespace PSEP {
  /*
   * SolutionProtocol is an enum class for one of the two arguments that may
   * be passed to the TSPSolver constructor
   * See pivplan.h for tuning the parameters of each of these protocols
   *
   * PURECUT - Attempts to solve the instance by adding only primal cutting
   *    planes
   * ABC - Runs PURECUT protocol until some specified termination condition,
   *    then embeds it in an Augment-Branch-Cut solver
   */
  enum class SolutionProtocol {
    PURECUT, ABC
      };

  /*
   * A very simple structure for creating a randomly generated euclidean TSP
   * instance.
   *
   * nodecount - number of cities. a value of zero shall be considered an empty
   *    or 'null' problem
   * gridsize - specifies the range of coordinates possible
   * seed - the random seed. if set to zero, problem will be generated 
   *    randomly using the current time as seed. Nonzero values will reliably
   *    generate the same problem.
   */
  struct RandProb {
  RandProb() : nodecount(0), gridsize(100), seed(0) {}
    int nodecount;
    int gridsize;
    int seed;
  };

  /*
   * The LP namespace here 
   * contains constants, enum classes, and stuctures related 
   * to the LP solver/solutions/relaxations. 
   */
  namespace LP {

    /* Tolerance by which two doubles are considered equal, or a number is
     * considered equal to zero
     */
    constexpr double EPSILON = 0.000001;

    /* The CPLEX default iteration limit for simplex optimizers */
    constexpr long long DEFAULT_ITLIM = 9223372036800000000;

    /*
     * Pricing is an enum class for methods by which pivot variables may be
     * selected in the primal simplex method. See CPLEX documentation for more
     * info
     *
     * Devex
     * SlackSteepest is steepest edge with slack initial norms
     * Steepest is true steepest edge
     */
    enum class Pricing {
      Devex, SlackSteepest, Steepest
    };

    /*
     * PivType is how we categorize LP solutions 
     *
     * Frac - fractional solution
     * Subtour - integral subtour
     * Tour - an augmented integral tour
     * FathomedTour - an integral tour that is dual feasible, i.e., optimal
     *    for the current LP relaxation. 
     *    With PURECUT solution method, fathomed is synonymous with optimality
     *    With ABC, this only means optimal at the current branching node
     */
    enum class PivType {
      Frac, Subtour, Tour, FathomedTour
    };

    /*
     * Prefs is a simple struct to store preferences for the LP solver.
     *
     * price_method is one of the Pricing types indicated above
     * dp_threshold - controls when simple domino parity inequality separation
     *    will be called. Negative values mean no simple DP separation, no 
     *    matter what. Nonnegative values mean simple DP separation will be
     *    called only after that many rounds of cutting plane generation, or
     *    if no other cuts are found
     */
    struct Prefs {
    Prefs() : price_method(Pricing::Devex), dp_threshold(-1) {}
    Prefs(Pricing _price, int _dp_threshold) :
      price_method(_price), dp_threshold(_dp_threshold) {}
      
      Pricing price_method;
      int dp_threshold;
    };
  }


  /*
   * Some timing functions: zeit for CPU time, real_zeit for wall clock time
   */
  double zeit (void);
  double real_zeit (void);

  //utility function to print every entry of a vector
  template<typename entry_t>
    void print_vec(std::vector<entry_t> const &vec){
    for(int i = 0; i < vec.size(); i++)
      std::cout << "Entry " << i << ": " << vec[i] << "\n";
  }

  //utility function to print nonzero entries of a numeric vector
  template<
    typename entry_t,
    typename =
    typename std::enable_if<std::is_arithmetic<entry_t>::value, entry_t>::type
    >
    void print_vec_nonzero(std::vector<entry_t> const &vec){
    for(int i = 0; i < vec.size(); i++)
      if(vec[i] != 0)
	std::cout << "Entry " << i << ": " << vec[i] << "\n";
  }
    
}

typedef std::pair<int, int> IntPair;
typedef std::unordered_map<IntPair, int, boost::hash<IntPair>> IntPairMap;


#endif
