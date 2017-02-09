/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**
 * @file
 * @brief Utility functions, macros, and structures.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_UTIL_H
#define CMR_UTIL_H

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <cmath>

/// The namespace for this project.
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


/// Numerical tolerances used in this project.
namespace Epsilon {


constexpr double Zero = 0.000001; //!< Numbers less than this treated as zero.

/// Cuts are not considered violated unless violated by at least this much.
constexpr double Cut = 0.0001;

/// Try another separation routine if pivot delta is less than this much.
constexpr double SepRound = 0.0001;

}

/// Utility functions/structures used miscellaneous places in the project.
namespace util {

/// Is a zero-one variable considered integral.
inline bool var_integral(double d)
{ return fabs(d) < Epsilon::Zero || fabs(d) > 1 - Epsilon::Zero; }

double zeit (void); //!< CPU time function.
double real_zeit (void); //!< Wall clock time function.

/// As per Herb Sutter, port of C++14's make_unique faculty. 
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
    /// Call operator for deleting the managed pointer.
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
template<typename numtype>
using c_array_ptr = std::unique_ptr<numtype, C_resource_deleter<numtype>>;

/// Class template for a square upper triangular matrix.
template<class T>
class SquareUT {
public:
    SquareUT() = default; //!< Default construct an empty matrix.

    /// Move constructor.
    SquareUT(SquareUT &&M) noexcept : data(std::move(M.data)) {}

    /// Move assignment.
    SquareUT &operator=(SquareUT &&M) noexcept
        { data = std::move(M.data); return *this; }
    

    /// Construct a matrix of dimension \p size.
    SquareUT(size_t size) : data(size)
        {
            for (auto i = 0; i < size; ++i)
                data[i] = std::vector<T>(size - i);
        }

    /// Construct a matrix of dimension \p size with \p val as all entries. 
    SquareUT(size_t size, T val) : data(size)
        {
            for (auto i = 0; i < size; ++i)
                data[i] = std::vector<T>(size - i, val);
        }

    /// Access entry by matrix subscripting. 
    T &operator()(size_t row, size_t column)
        {
            return data[row][column - row];
        }

private:
    /// The entries of the matrix, stored as a ragged vector of vectors.
    std::vector<std::vector<T>> data;
};

}

constexpr int IntMax = std::numeric_limits<int>::max();
constexpr int IntMin = std::numeric_limits<int>::min();

constexpr double DoubleMax = std::numeric_limits<double>::max();
constexpr double DoubleMin = std::numeric_limits<double>::min();

typedef std::pair<int, int> IntPair;

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

/// Simple base class for storing edge of a graph as a sorted pair of nodes.
struct EndPts {
    EndPts() = default;
    EndPts(int e0, int e1) : end{{e0, e1}}
        { if (end[0] > end[1]) std::swap(end[0], end[1]); }

    std::array<int, 2> end;
};
    
}




#endif
