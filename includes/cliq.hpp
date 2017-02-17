/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Reference counted storage of Cliques and Tooth objects.
 * The structures in this file are meant to be managed by a HyperGraph
 * and ExternalCuts.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_CLIQ_H
#define CMR_CLIQ_H

#include "cut_structs.hpp"
#include "util.hpp"

#include <array>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace CMR {
namespace Sep {

/** Class for storing segment lists representing edges of a hypergraph.
 * A Clique stores a subset of vertices as a list of CMR::Segment objects,
 * where the start and endpoints indicate a range of nodes from a tour.
 * Thus, a Clique is meaningless without a tour from which to be
 * derefrenced.
 */
class Clique {
public:
    Clique() = default; //!< Default construct an empty Clique.

    /// Construct a Clique from a Concorde clique.
    Clique(const CCtsp_lpclique &cc_cliq,
           const std::vector<int> &saved_tour,
           const std::vector<int> &saved_perm,
           const std::vector<int> &current_tour);

    /// Construct a Clique from start and end indices.
    Clique(int start, int end,
           const std::vector<int> &saved_tour,
           const std::vector<int> &saved_perm,
           const std::vector<int> & current_tour);

    /// Construct a Clique from a list of literal nodes.
    Clique(std::vector<int> &nodes, const std::vector<int> &perm);

    using Ptr = std::shared_ptr<Clique>; //!< shared_ptr alias declaration.

    /// How many segments are used to represent the Clique.
    int seg_count() const { return seglist.size(); }

    /// A constant reference to the list of segments in the Clique.
    const std::vector<Segment> &seg_list() const { return seglist; }

    /// A list of literal nodes represented by the Clique.
    std::vector<int> node_list(const std::vector<int> &saved_tour) const;

    bool operator==(const Clique &rhs) const
        { return seglist == rhs.seglist; } //!< Equality operator.

    /// Returns true iff the Clique contains \p index.
    bool contains(const int index) const {
        for (const Segment &seg : seglist)
            if (seg.contains(index))
                return true;
        return false;
    }

private:
    /// A vector of start and endpoints of tour intervals stored as Segment.
    std::vector<Segment> seglist;
};

/** Vertex set structure used in tooth inequalities for domino parity cuts.
 * This class holds a tooth in the sense used by Fleischer et al. (2006)
 * meaning two disjoint, nonempty vertex subsets whose union is not the vertex
 * set of the graph. A tooth inequality is obtained by summing the SECs on the
 * two sets.
 * @remark This underlying data in this class is sufficiently general to hold
 * general teeth (i.e., for general domino parity inequalities), but for now
 * there is only a constructor to implement simple teeth for simple domino
 * parity inequalities, i.e., teeth where one of the sets is a singleton.
 */
class Tooth {
public:
    Tooth() = default; //!< Default construct an empty tooth.

    /// Construct a Tooth from a SimpleTooth.
    Tooth(const SimpleTooth &T,
          const std::vector<int> &saved_tour,
          const std::vector<int> &saved_perm,
          const std::vector<int> &current_tour);

    /// Constant reference to the defining sets.
    const std::array<Clique, 2> &set_pair() const { return sets; }

    using Ptr = std::shared_ptr<Tooth>; //!< Pointer alias declaration.

    bool operator==(const Tooth &rhs) const
        { return sets == rhs.sets; } //!< Equality operator.

private:
    /// Teeth are represented as a pair of Cliques for the handle and body.
    std::array<Clique, 2> sets;
};


}
}

namespace std {

/// Partial specialization of std::hash taken from CCtsp_hashclique.
template<>
struct hash<CMR::Sep::Clique> {
    /// Call operator for hashing a Clique.
    size_t operator()(const CMR::Sep::Clique &clq) const
        {
            size_t val = 0;

            for (const CMR::Segment &seg : clq.seg_list())
                val = (val * 65537) + (seg.start * 4099) + seg.end;

            return val;
        }
};

/// Partial specialization of std::hash from CCtsp_hashdomino.
template<>
struct hash<CMR::Sep::Tooth> {
    /// Call operator for hashing a Tooth.
    size_t operator()(const CMR::Sep::Tooth &T) const
        {
            size_t val = 0;

            for (const CMR::Sep::Clique &clq : T.set_pair())
                for (const CMR::Segment &seg : clq.seg_list())
                    val = (val * 65537) + (seg.start * 4099) + seg.end;

            return val;
        }
};

}

namespace CMR {
namespace Sep {

/** Storage of a repository of Cliques, for use in building a HyperGraph.
 * This class is responsible for dispensing and deleting references to Clique
 * objects, for use in the edges of a HyperGraph. The interface allows
 * a Clique to be passed to the CliqueBank. A pointer to the Clique will be
 * added to the bank if one does not exist already, or else the use count
 * of the pointer will be incremented. This pointer is then returned for other
 * use. The Cliques contained
 * therein are meaningless without reference to a fixed tour and perm vector,
 * which shall be used to construct Clique objects and turn them back into
 * lists of nodes or sparse cut rows.
 */
class CliqueBank {
public:
    /// Construct a CliqueBank to be dereferenced by \p tour and \p perm.
    CliqueBank(const std::vector<int> &tour, const std::vector<int> &perm);

    /// Add a Clique to the bank, and get a reference to it.
    Clique::Ptr add_clique(const Clique &clq);

    /// Construct a Clique in place, add it, and get a reference to it.
    Clique::Ptr add_clique(const CCtsp_lpclique &cc_cliq,
                                const std::vector<int> &tour);

    /// Construct/add/get a reference to the Clique from endpoints.
    Clique::Ptr add_clique(int start, int end,
                                const std::vector<int> &tour);

    /// Construct/add/get a reference to a Clique from a node list.
    Clique::Ptr add_clique(std::vector<int> &nodes);

    /// Put the pointed Clique in this bank.
    void steal_clique(Clique::Ptr &clq_ptr, CliqueBank &from_bank);

    /// Decrement the reference count of a Clique, possibly removing it.
    void del_clique(Clique::Ptr &clq_ptr);

    int size() const
        { return bank.size(); } //!< The number of unique Cliques in the bank.

    /// Alias declaration for Clique hash table.
    using CliqueHash = std::unordered_map<Clique, Clique::Ptr>;

    using Itr = CliqueHash::iterator;
    using ConstItr = CliqueHash::const_iterator;

    Itr begin() { return bank.begin(); }
    Itr end() { return bank.end(); }

    ConstItr begin() const { return bank.begin(); }
    ConstItr end() const { return bank.end(); }

    const std::vector<int> &ref_tour() const { return saved_tour; }
    const std::vector<int> &ref_perm() const { return saved_perm; }


private:
    const std::vector<int> saved_tour; //!< Saved tour for dereferencing.
    const std::vector<int> saved_perm; //!< Permutation vector for saved_tour.
    CliqueHash bank; //!< Hash table of Clique and Clique::Ptr.
};


/** Storage of a repository of Teeth, for use in building a HyperGraph.
 * This class is like CliqueBank, but instead it dispenses and deletes
 * references to Tooth objects instead.
 */
class ToothBank {
public:
    /// Construct a ToothBank to be dereferenced by \p tour and \p perm.
    ToothBank(const std::vector<int> &tour, const std::vector<int> &perm);

    /// Construct a ToothBank with reference vectors matching a CliqueBank.
    ToothBank(const CliqueBank &cbank);

    /// Add a Tooth to the bank and get a reference to it.
    Tooth::Ptr add_tooth(const SimpleTooth &T,
                         const std::vector<int> &tour);

    /// Decrement the reference count of a Tooth, possibly deleting it.
    void del_tooth(Tooth::Ptr &T_ptr);

    int size() const
        { return bank.size(); } //!< Number of teeth in the bank.

    /// Alias declaration for Tooth hash table.
    using ToothHash = std::unordered_map<Tooth, Tooth::Ptr>;

    using Itr = ToothHash::iterator; //!< Iterator alias.
    using ConstItr = ToothHash::const_iterator; //!< Const iterator alias.

    Itr begin()
        { return bank.begin(); } //!< Begin iterator.

    ConstItr begin() const
        { return bank.begin(); } //!< Const begin iterator.

    Itr end()
        { return bank.end(); } //!< Past the end iterator.

    ConstItr end() const
        { return bank.end(); } //!< Const past the end iterator.

    const std::vector<int> &ref_tour() const
        { return saved_tour; } //!< Const ref to saved tour for dereferencing.

    const std::vector<int> &ref_perm() const
        { return saved_perm; } //!< Const ref to saved perm for dereferencing.

private:
    const std::vector<int> saved_tour; //!< Saved tour for dereferencing.
    const std::vector<int> saved_perm; //!< Permutation vec for saved_tour.
    ToothHash bank; //!< Hash table of Tooth to Tooth::Ptr.
};

}
}

#endif
