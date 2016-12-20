/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief UTILITY FUNCTIONS AND MACROS/ENUM CLASSES
 *
 * This header contains namespace and enum classes for labels/constants that
 * are used elsewhere in the code, as well as some very simple structures
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _CMR_UTIL_H
#define _CMR_UTIL_H

#include <utility>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <limits>
#include <string>
#include <memory>

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

/** Function for converting PivType to a string for output. */
std::string piv_string(PivType piv);

}


/** Namespace for numerical tolerances used in this project. */
namespace Epsilon {

/** Numbers less than this are treated as zero. */
constexpr double Zero = 0.000001;

/** Cuts are not considered violated unless violated by at least this much. */
constexpr double Cut = 0.0001;

}

/** CPU time function. */
double zeit (void);

/** Wall-clock time function. */
double real_zeit (void);

/** Exception-safe creation of a unique_ptr, from Herb Sutter.
 * This is a c++14 feature which is ported to c++11 under the CMR namespace,
 * since it exists as std::make_unique in c++14 and onwards.
 */
template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args )
{
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}

/** Class template for deleting resources allocated by C functions.
 * In calls to Concorde, memory allocation takes place through malloc and 
 * free, and it can be dangerous to delete memory that was malloc'd. Thus,
 * A Concorde resource managed by a unique_ptr shall use this deleter to 
 * free the memory associated with the resource.
 */
template <typename T>
struct C_resource_deleter {
    void operator()(T* ptr) const {
        if (ptr) free(ptr);
        ptr = nullptr;
    }
};


/** Alias declaration for unique_ptr to C array.
 * This specialization of unique_ptr takes ownership of a C-style array which
 * was dynamically allocated by a C routine, freeing its memory 
 * appropriately when necessary.
 */
using c_array_ptr = std::unique_ptr<int, C_resource_deleter<int>>;

/** Alias template for C structs with free/init functions.
 * Suggested usage is as follows. In Concorde there are many structures of the
 * form `struct CCtsp_struct` which are used as `CCtsp_struct *S`
 * which is malloc'd, initialized by `CCtsp_init_struct(S)`, and then the
 * associated memory is freed by a call to `CCtsp_free_struct(S)`. Hence
 * declare `*S` and call `CCtsp_init_struct(S)`, and then manage the object
 * by c_struct_ptr(S, &CCts_free_struct) to create a unique_ptr which will
 * automatically free S using the appropriate Concorde function to release
 * the memory. 
 */
template <typename T>
using c_struct_ptr = std::unique_ptr<T, void(*)(T*)>;

constexpr int IntMax = std::numeric_limits<int>::max();
constexpr int IntMin = std::numeric_limits<int>::min();

constexpr double DoubleMax = std::numeric_limits<double>::max();
constexpr double DoubleMin = std::numeric_limits<double>::min();

typedef std::pair<int, int> IntPair;

struct edge_hash {
    std::size_t operator()(const IntPair& p) const
        {
            std::hash<int> hasher;
            std::size_t seed = 0;
            seed ^= hasher(p.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hasher(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed; 
        }
};

typedef std::unordered_map<IntPair, int, edge_hash> IntPairMap;

/** Simple utility struct for storing an interval of nodes.
 * A Segment is defined in terms of some list of nodes, usually a tsp tour.
 * If \p tour is the nodelist, then a Segment `S` defined relative to tour
 * represents the nodes `tour[S.start], ..., tour[S.end]`. Thus a Segment
 * is meaningless without a tour from which to be dereferenced. 
 * @remark Representation of nodes as a segment of a tour is a common theme
 * in tsp computation, and for primal cutting plane tsp computation in 
 * particular. This class is open for inheritance to define subtour cuts 
 * associated to segments, bodies of simple teeth, and simple teeth themselves
 * where the notion of size and containment is identical.
 */
struct Segment {
    /** Default construct a Segment. */
    Segment() = default;

    /** Construct a Segment with specified start and end point. */
    Segment(int lo, int hi) : start(lo), end(hi) {}

    /** Size of the Segment.
     * This is the number of the nodes in the closed interval from 
     * \p start to \p end.
     */
    int size() const { return end - start + 1; }

    /** Does the Segment contain a certain vertex.
     * Returns true iff \p vx lies in the interval specified by the Segment.
     */
    bool contains(int vx) const { return start <= vx && vx <= end; }

    /** Is one Segment a subset of the other.
     * Returns true iff the interval defined by this Segment is a subset of
     * that defined by \p seg.
     */
    bool subset_of(const Segment &seg) const
    { return seg.contains(start) && seg.contains(end); }

    /** Compare segments by size and then start point. */
    bool operator>(const Segment &rhs) const
    {
        return std::make_tuple(size(), start, end) >
        std::make_tuple(rhs.size(), rhs.start, rhs.end);
    }

    /** Equality operator. */
    bool operator==(const Segment &rhs) const
    { return start == rhs.start && end == rhs.end; }

    int start; /**< The start index of the Segment. */
    int end; /**< The end index of the Segment. */
};
    
}




#endif
