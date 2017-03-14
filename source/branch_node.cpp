#include "branch_node.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

using std::runtime_error;
using std::logic_error;
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
        result << "[" << B.ends << " = " << B.direction << "]";
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
        os << "depth " << B.depth << ", tourlen "
           << static_cast<int>(B.tourlen)
           << ", parent "
           << bnode_brief(*(B.parent));
    }
    return os;
}


}
}
