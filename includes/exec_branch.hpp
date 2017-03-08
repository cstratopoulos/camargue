/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Branching execution.
 */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_EXEC_BRANCH_H
#define CMR_EXEC_BRANCH_H

#include "active_tour.hpp"
#include "core_lp.hpp"
#include "datagroups.hpp"
#include "cliq.hpp"
#include "lp_util.hpp"
#include "branch_util.hpp"

#include <array>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

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
               Sep::Clique::Ptr tour_clq_, double tourlen_);

    BranchNode(BranchNode &&B) noexcept;
    BranchNode &operator=(BranchNode &&B) noexcept;

    /// Alias declaration for returning two split child problems.
    using Split = std::array<BranchNode, 2>;

    EndPts ends; //!< The endpoints of the branch edge.
    Dir direction; //!< Down branch or up branch.

    Status stat; //!< The type of processing required by the node.

    const BranchNode *parent;
    int depth; //!< Search tree depth of this node.

    Sep::Clique::Ptr tour_clq; //!< The tour to be instated at this node.
    double tourlen; //!< The length of the tour in tour_clq.

    /// A starting basis for if Status is NeedsPrice or NeedsRecover.
    LP::Basis::Ptr price_basis;

    /// Is this the root problem.
    bool is_root() const { return parent == nullptr; }

    /// Has the problem been processed.
    bool visited() const
        { return stat == Status::Pruned || stat == Status::Done; }
};

/// Turn a 0-1 value \p entry into a BranchNode::Dir.
BranchNode::Dir dir_from_int(int entry);

/// A concise "label" for the BranchNode \p B.
std::string bnode_brief(const BranchNode &B);

std::ostream &operator<<(std::ostream &os, const BranchNode::Dir &dir);
std::ostream &operator<<(std::ostream &os, const BranchNode &B);

using BranchHistory = std::list<BranchNode>;
using SplitIter = std::array<BranchHistory::iterator, 2>;


class Executor {
public:
    /// Construct an Executor using data from an existing solution process.
    Executor(const Data::Instance &inst, const LP::ActiveTour &activetour,
             const Data::BestGroup &bestdata,
             const Graph::CoreGraph &coregraph, LP::CoreLP &core);

    /// Alias declaration for EndPts and branching direction.
    using EndsDir = std::pair<EndPts, BranchNode::Dir>;

    ScoreTuple branch_edge(); //!< Get the next edge to branch on.

    /// Create the children nodes of \p parent for branching on \p branch_edge.
    BranchNode::Split split_problem(ScoreTuple &branch_tuple,
                                    BranchNode &parent);

    /// Get a tour satisfying the branching in \p edge_stats.
    void branch_tour(const std::vector<EndsDir> &edge_stats,
                     const std::vector<int> &start_tour_nodes,
                     bool &found_tour, bool &feas,
                     std::vector<int> &tour, double &tour_val);

    /// Compress \p tour into a Clique reference.
    Sep::Clique::Ptr compress_tour(const std::vector<int> &tour);

    /// Expand a compressed representation computed by compress_tour.
    std::vector<int> expand_tour(const Sep::Clique::Ptr &tour_clique);

    /// Clamp a variable as indicated by \p current_node.
    void clamp(const BranchNode &current_node);

    /// Undo the clamp done on \p current_node.
    void unclamp(const BranchNode &current_node);

private:
    const Data::Instance &instance;
    const LP::ActiveTour &active_tour;
    const Data::BestGroup &best_data;
    const Graph::CoreGraph &core_graph;

    LP::CoreLP &core_lp;

    /// Used by branch_tour to track the edges in \p edge_stats.
    /// An instance.node_count() dimensional upper triangular matrix with
    /// `inds_table(i, j)` being 'A' if the edge is to be avoided, 'W' if
    /// the edge is wanted (fixed to one), and '\0' otherwise.
    util::SquareUT<char> inds_table;

    /// Used by branch_tour to track fixed edges.
    /// An `instance.node_count()` length array initialized to all zero, and
    /// incremented for the ends of every BranchNode::Dir::Up edge in
    /// edge_stats. No tour can exist if any entry is greater than two.
    std::vector<int> fix_degrees;

    /// Used by branch_tour to track avoided edges.
    /// Like fix_degrees, but for track BranchNode::Dir::Down edges in
    /// edge_stats. Initialized to all `instance.node_count() - 1`, and
    /// decremented. No tour can exist if an entry is less than two.
    std::vector<int> avail_degrees;

    /// Used to compress the tours returned by branch_tour.
    /// This CliqueBank shall represent Cliques in terms of whatever was the
    /// best_data.best_tour_nodes when this Executor was constructed.
    Sep::CliqueBank tour_cliques;
};

}
}

#endif
