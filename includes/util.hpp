/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief UTILITY FUNCTIONS AND MACROS/ENUM CLASSES
 *
 * This header contains namespace and enum classes for labels/constants that
 * are used elsewhere in the code, as well as some very simple structures
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef __CMR_UTIL_H
#define __CMR_UTIL_H

#include <utility>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <type_traits>
#include <limits>
#include <string>

#include <boost/functional/hash.hpp>

/** Macro for printing error message and going to exit label. 
 * Prints the message \a message to stderr and goes to CLEANUP. Thus CLEANUP
 * must have been defined as a labelled section in function scope.
 */
#define CMR_GOTO_CLEANUP(message) { std::cerr << message; goto CLEANUP; }

/** Conditional version of #CMR_GOTO_CLEANUP(message).
 * Only prints the error \a message if \a rval is nonzero.
 */
#define CMR_CHECK_RVAL(rval, message) {	\
    if ((rval)) {				\
      std::cerr << message;			\
      goto CLEANUP;				\
    }						\
  }

/** Version of #CMR_GOTO_CLEANUP(message) for setting return code.
 * This macro sets rval to 1, prints \a message to stderr, and goes to 
 * CLEANUP. Its intended use is in mixing try/catch blocks with retcodes, and
 * in function calls where subtasks have their own retcodes. For example,
 * `try { task } catch (...) { CMR_SET_GOTO(rval, "task threw exception"); }`
 * or `subval = sub_task(); if(subval) CMR_SET_GOTO(rval, "sub_task error")`
 */
#define CMR_SET_GOTO(rval, message) {		\
    rval = 1;					\
    std::cerr << message;			\
    goto CLEANUP; }

/** The namespace for this project. */
namespace CMR {

/** Choices for the execution behavior of PureCut. */
enum class SolutionProtocol {

  /** Pure cutting plane solution
   * Indicates attempt to solve solely by pivoting and adding primal cutting
   * planes, continuing until optimality or failure. 
   */
  PURECUT,

  /** Augment-Branch-Cut solution
   * Indicates that `PureCut.solve` will run until some specified threshold of
   * time or cuts, and then will be embedded in an Augment-and-Branch-and-Cut
   * solution process if optimality has not been attained.
   */
  ABC
};

/** A POD struct for creating a randomly generated euclidean TSP instance. */
struct RandProb {
  RandProb() : nodecount(0), gridsize(100), seed(0) {}
  int nodecount; /**< Number of cities in the problem. */
  int gridsize; /**< Size of grid from which coordinates will be drawn. */
  int seed; /**<Random seed for generation of points. */
};

/** A POD struct for preferences related to tour file output.
 * This structure stores preferences about the types of files (if any) that
 * will be saved to disk during the solution process. Existing files will
 * be overwritten; it is the responsibility of the calling routine to make
 * sure that they are hopefully only overwritten by better tours than before!
 * @todo Add a make GIF type parameter. 
 */
struct OutPrefs {
  OutPrefs() : save_tour(true),
	       dump_xy(false),
	       save_tour_edges(false),
	       probname() {}
  OutPrefs(const bool _save_tour, const bool _dump_xy,
	   const bool _save_tour_edges) :
    save_tour(_save_tour),
    dump_xy(_dump_xy),
    save_tour_edges(_save_tour_edges),
    probname() {}
  
  bool save_tour, /**< Save tour nodes to `probname.sol`. */
    dump_xy, /**< If possible, dump xy-coords to `probname.xy` */
    save_tour_edges; /**< Save tour edges to `probname_tour.x`. */
  std::string probname; /**< The name of the problem. */
};

/** Namespace for classes, constants, and enums related to %LP relaxations. */
namespace LP {

/** The CPLEX default iteration limit for simplex optimizers */
constexpr long long DEFAULT_ITLIM = 9223372036800000000;

/** Enum class for primal pricing protocols. 
 * Pricing is an enum class for methods by which pivot variables may be
 * selected in the primal simplex method. See CPLEX documentation for more
 * info. These appear in order of least to most computationally intensive and,
 * theoretically, least to most effective at identifying non-degenerate pivots.
 */
enum class Pricing {
  Devex, /**< Devex pricing. */
  SlackSteepest, /**< Steepest edge with slack initial norms. */
  Steepest /**< True steepest edge. */
};

/** Enum class for categorizing lp solutions. */
enum class PivType {
  Frac, /**< Fractional solution. */
  Subtour, /**< Integral subtour. */
  Tour, /**< A new or augmented tour. */
  FathomedTour /**< A Tour with a dual feasible basis in the current lp. */
};

std::string piv_string(PivType piv);

/*
 * Prefs is a simple struct to store preferences for the lp solver.
 *
 * price_method is one of the Pricing types indicated above
 * dp_threshold - controls when simple domino parity inequality separation
 *    will be called. Negative values mean no simple DP separation, no 
 *    matter what. Nonnegative values mean simple DP separation will be
 *    called only after that many rounds of cutting plane generation, or
 *    if no other cuts are found
 */

/** POD struct for lp solver preferences. */
struct Prefs {
  Prefs();
  Prefs(Pricing _price, int _dp_threshold, int max_round, int q_max);
      
  Pricing price_method; /**< Which primal pricing criterion to use. */
      
  int dp_threshold; /**< Wait this many rounds before calling simple DP sep. */
  int max_per_round; /**< Add at most this many cuts of each type per round. */
  int q_max_size; /**< Keep a pool of this many blossoms. */
};
}


/** Namespace for numerical tolerances used in this project. */
namespace Epsilon {

/** Numbers less than this are treated as zero. */
constexpr double Zero = 0.000001;

/** Cuts are not considered violated unless violated by at least this much. */
constexpr double Cut = 0.0001;

/** The connected component tolerance for Grotschel Holland fast blossoms. */
constexpr double GH = 0.3;

}

/** CPU time function. */
double zeit (void);

/** Wall-clock time function. */
double real_zeit (void);

/** Prints every entry of a vector */
template<typename entry_t>
void print_vec(std::vector<entry_t> const &vec){
  for(int i = 0; i < vec.size(); i++)
    std::cout << "Entry " << i << ": " << vec[i] << "\n";
}

/** Prints nonzero entries of a vector. */
template<typename entry_t >
void print_vec_nonzero(std::vector<entry_t> const &vec){
  for(int i = 0; i < vec.size(); i++)
    if(vec[i] != 0)
      std::cout << "Entry " << i << ": " << vec[i] << "\n";
}

/** Exception-safe creation of a unique_ptr, from Herb Sutter.
 * This is a c++14 feature which is ported to c++11 under the CMR namespace,
 * since it exists as std::make_unique in c++14 and onwards.
 */
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

template <typename T>
struct C_resource_deleter {
    void operator()(T* ptr) const {
        if(ptr) free(ptr);
        ptr = nullptr;
    }
};

typedef std::unique_ptr<int, C_resource_deleter<int>> c_array_ptr;

constexpr int IntMax = std::numeric_limits<int>::max();
constexpr int IntMin = std::numeric_limits<int>::min();

constexpr double DoubleMax = std::numeric_limits<double>::max();
constexpr double DoubleMin = std::numeric_limits<double>::min();
    
}

typedef std::pair<int, int> IntPair;
typedef std::unordered_map<IntPair, int, boost::hash<IntPair>> IntPairMap;


#endif
