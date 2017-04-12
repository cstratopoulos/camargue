#include "branch_node.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <tuple>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

using std::runtime_error;
using std::exception;

namespace CMR {
namespace ABC {

constexpr int IntMax = std::numeric_limits<int>::max();
constexpr double DoubleMax = std::numeric_limits<double>::max();

BranchNode::BranchNode() : stat(Status::NeedsCut),
                           parent(nullptr), depth(0),
                           tourlen(IntMax),
                           estimate(DoubleMax) {}

BranchNode::BranchNode(EndPts ends_, Dir direction_,
                       const BranchNode &parent_,
                       double tourlen_, double estimate_)
    : ends(ends_), direction(direction_), stat(Status::NeedsCut),
      parent(&parent_), depth(1 + parent_.depth),
      tourlen(tourlen_), estimate(estimate_) {}

BranchNode::BranchNode(BranchNode &&B) noexcept
    : ends(std::move(B.ends)),
      direction(B.direction),
      stat(B.stat),
      parent(B.parent),
      depth(B.depth),
      tourlen(B.tourlen),
      price_basis(std::move(B.price_basis)),
      estimate(B.estimate)
{
    B.stat = Status::Done;

    B.parent = nullptr;
    B.depth = 0;

    B.tourlen = IntMax;
    B.estimate = DoubleMax;
}

BranchNode &BranchNode::operator=(BranchNode &&B) noexcept
{
    if (&B == this)
        return *this;

    ends = std::move(B.ends);
    direction = B.direction;

    stat = B.stat;


    parent = B.parent;
    depth = B.depth;

    tourlen = B.tourlen;

    estimate = B.estimate;
    price_basis = std::move(B.price_basis);

    B.stat = Status::Done;

    B.parent = nullptr;
    B.depth = 0;

    B.tourlen = IntMax;
    B.estimate = DoubleMax;

    return *this;
}

/**
 * Comparator where \p A is worse than \p B if \p A is estimated to have a
 * longer best tour. Ties are broken by depth and then strong branch estimate,
 * in that order.
 * @remark the depth of \p A and \p B are reversed `tie`d because we prefer
 * greater values of depth but lesser values of tour and estimate.
 */
bool BranchNode::tour_worse(const BranchHistory::iterator &A,
                            const BranchHistory::iterator &B)
{
    return (std::tie(A->tourlen, B->depth, A->estimate) >
            std::tie(B->tourlen, A->depth, B->estimate));
}

/**
 * Comparator where \p A is worse than \p B if \p A is estimated to have a
 * higher objective value lower bound. Ties are broken by depth and estimated
 * tour length, in that order.
 * @remark the depth of \p A and \p B are reversed `tie`d because we prefer
 * greater values of depth but lesser values of tour and estimate.
 */
bool BranchNode::bound_worse(const BranchHistory::iterator &A,
                             const BranchHistory::iterator &B)
{
    return (std::tie(A->estimate, B->depth, A->tourlen) >
            std::tie(B->estimate, A->depth, B->tourlen));
}

BranchNode::Dir dir_from_int(int entry)
{
    switch (entry) {
    case 0:
        return BranchNode::Dir::Down;
    case 1:
        return BranchNode::Dir::Up;
    default:
        throw runtime_error("Tried to get BranchNode::Dir w non-binary int");
    }
}

std::string bnode_brief(const BranchNode &B)
{
    std::stringstream result;

    if (B.is_root())
        result << "[root]";
    else {
        result << "[" << B.ends << " = " << B.direction << "], lvl "
               << B.depth;
    }

    return result.str();
}

ostream &operator<<(ostream &os, const BranchNode::Dir &dir)
{
    os << static_cast<int>(dir);
    return os;
}
ostream &operator<<(ostream &os, const BranchNode &B)
{
    os << bnode_brief(B) << ", ";

    if (B.is_root()) {
        os << "depth 0";
    } else {
        os << "tour "
           << static_cast<int>(B.tourlen)
           << ", est " << std::setprecision(2) << B.estimate
           << std::setprecision(6) << std::fixed
           << ", parent "
           << bnode_brief(*(B.parent));
    }
    return os;
}


}
}
