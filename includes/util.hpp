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
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <cmath>

/// The namespace for this project.
namespace CMR {

/// Preferences related to tour file output and verbosity.
struct OutPrefs {
    std::string probname; //!< The instance name.

    /**@name Tour/instance file output.
     * Each parameter if true determines frequency/style of file output.
     * @see util::write_tour_nodes for method and style of tour node output.
     * @see util::write_tour_edges for tour edge output.
     */
    ///@{

    bool save_tour = true; //!< Save the latest, best tour to `probname.sol`.
    bool dump_xy = false; //!< If possible, dump x-y coords to `probname.xy`
    bool save_tour_edges = false; //!< Save save_tour edges to `probname_tour.x`

    /// Record each augmenting tour as specified by save_tour/save_tour_edges.
    /// If true, then for each file output will have an infix indicating what
    /// number augmentation this is. For example if save_tour_edges is true
    /// and this value is true, a sequence a files `probname_tour.0.x`,
    /// `probname_tour.1.x`, etc. will be created.
    bool gif_tour = false;

    ///@}

    bool verbose = false; //!< Verbose output from classes and subroutines.

    /// A jittering progress bar for pivot values.
    bool prog_bar = false;

    /// Detailed timer/profiling of the code.
    bool detailed_stats = false;
};


/// Numerical tolerances.
namespace Epsilon {


constexpr double Zero = 0.000001; //!< Numbers less than this treated as zero.
constexpr double MinCut = 0.0001; //!< Tolerance for min cut computations.

/// Cuts are not considered violated unless by at least this much.
constexpr double CutViol = 0.001;

constexpr double DualDust = 0.001; //!< Small dual values.

/// A round of cuts is a failure if the pivot deltas sum to less than this.
constexpr double TotalDelta = 0.01;

constexpr double PHratio = 0.1; //!< A small value of the Padberg-Hong metric.

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

/// Adaptation of make_unique to reset a unique_ptr to an object.
/// For classes with reference members, we cannot easily define move
/// constructors. This function provides a surrogate move constructor by
/// resetting a unique_ptr to such an object.
template<typename T, typename ...Args>
void ptr_reset(std::unique_ptr<T> &target, Args&& ...args)
{
    target.reset(new T(std::forward<Args>(args)...));
}

/// Standardization of erase-remove idiom.
/// @param[in, out] vec the vector to modify.
/// @param[in] pred the predicate to apply.
/// @post deletes all elements of \p vec for which \p pred returns true.
template<typename ElemType, typename PredType>
void erase_remove(std::vector<ElemType> &vec, PredType pred)
{
    vec.erase(std::remove_if(vec.begin(), vec.end(), pred), vec.end());
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

    /// Fill all entries of this matrix with \p val.
    void fill(const T &val)
        {
            for (std::vector<T> &vec : data)
                std::fill(vec.begin(), vec.end(), val);
        }

private:
    /// The entries of the matrix, stored as a ragged vector of vectors.
    std::vector<std::vector<T>> data;
};

}

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
    EndPts()  : end{{-1, -1}} {}
    EndPts(int e0, int e1) : end{{e0, e1}}
        { if (end[0] > end[1]) std::swap(end[0], end[1]); }

    std::array<int, 2> end;
};

inline bool operator==(EndPts e1, EndPts e2)
{
    return e1.end == e2.end;
}

inline std::ostream &operator<<(std::ostream &os, const EndPts &e)
{
    os << "(" << e.end[0] << ", " << e.end[1] << ")";
    return os;
}

}




#endif
