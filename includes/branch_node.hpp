/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief ABC search nodes.
 * This file contains the structure representing branching subproblems to be
 * examined, and various comparator functions for implementing node selection
 * rules. There are also some I/O utility functions.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_NODE_H
#define CMR_BRANCH_NODE_H

#include "lp_util.hpp"
#include "util.hpp"

#include <array>
#include <iostream>
#include <list>
#include <utility>


namespace CMR {
namespace ABC {

struct BranchNode {
    enum Dir : int {
        Down = 0,
        Up = 1
    };

    enum class Status {
        NeedsCut,
        NeedsBranch,
        NeedsPrice,
        NeedsRecover,
        Pruned,
        Done,
    };

    BranchNode(); //!< Construct a root node.

    /// Construct a child node.
    BranchNode(EndPts ends_, Dir direction_, const BranchNode &parent_,
               double tourlen_, double estimate_);

    BranchNode(BranchNode &&B) noexcept;
    BranchNode &operator=(BranchNode &&B) noexcept;

    /// Alias declaration for returning two split child problems.
    using Split = std::array<BranchNode, 2>;

    EndPts ends; //!< The endpoints of the branch edge.
    Dir direction; //!< Down branch or up branch.

    Status stat; //!< The type of processing required by the node.

    const BranchNode *parent;
    int depth; //!< Search tree depth of this node.

    double tourlen; //!< Estimated best tour length for this node.

    /// A starting basis for if Status is NeedsPrice or NeedsRecover.
    LP::Basis::Ptr price_basis;
    double estimate; //!< The objective value estimate from edge selection.

    /// Is this the root problem.
    bool is_root() const { return parent == nullptr; }

    /// Has the problem been processed.
    bool visited() const
        { return stat == Status::Pruned || stat == Status::Done; }

    /**@name Node selection comparators.
     * These functions can be used to rank BranchNode objects so as to
     * implement a node selection rule.
     * @see abc_nodesel.hpp for the rules currently implemented.
     * @warning Other than depth-first search, node selection rules will be
     * implemented using std::make_heap, std::push_heap, std::pop_heap, etc.,
     * from the C++ STL. Such heaps are _max_ heaps by default, and node
     * selection rules usually demand _min_ heaps (e.g., shortest tour,
     * lowest strong branch bound). Thus, these functions must implement a
     * comparator where preferred nodes are _greater_. This should be
     * implemented as a lexicographic ordering with std::tie and/or
     *  std::make_tuple, with unvisited nodes always comparing greater than
     * visited ones, and then the actual ranking criteria. For example
     * if \f$ F \f$ is a comparator and \f$ A, B \f$ are two unvisited nodes
     *  where \f$ A \f$ has a shorter tour than \f$ B \f$, then
     * \f$ F(A, B) \f$ should return false.
     */
    ///@{

    /// Returns true if \p A has a longer tour than \p B.
    static bool tour_compare(const BranchNode &A, const BranchNode &B);

    /// Returns true if \p A has a higher strong branch estimate than \p B.
    static bool bound_compare(const BranchNode &A, const BranchNode &B);

    ///@}
};

/// Turn a 0-1 value \p entry into a BranchNode::Dir.
BranchNode::Dir dir_from_int(int entry);

/// A concise "label" for the BranchNode \p B.
std::string bnode_brief(const BranchNode &B);

std::ostream &operator<<(std::ostream &os, const BranchNode::Dir &dir);
std::ostream &operator<<(std::ostream &os, const BranchNode &B);

using BranchHistory = std::list<BranchNode>;

/// Alias declaration for EndPts and branching direction.
using EndsDir = std::pair<EndPts, BranchNode::Dir>;




}
}

#endif
