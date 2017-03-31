/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Working with branch tours.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_TOUR_H
#define CMR_BRANCH_TOUR_H

#include "branch_node.hpp"
#include "datagroups.hpp"
#include "graph.hpp"
#include "edgehash.hpp"

#include <utility>
#include <vector>


namespace CMR {
namespace ABC {

/// Compute branch tours for estimation and instatement, managing their edges.
class BranchTourFind {
public:
    BranchTourFind(const Data::Instance &inst,
                   Data::BestGroup &bestdata,
                   Graph::CoreGraph &coregraph);

    /// Compute a list of common constraints for splitting on \p parent.
    std::vector<EndsDir> common_constraints(const BranchNode &parent,
                                            const EndPts &branch_edge);

    /// Compute a list of branch constraints for \p N and all its ancestors.
    std::vector<EndsDir> branch_constraints(const BranchNode &N);

    /// Is the \p tour compliant with \p constraints.
    bool tour_compliant(const std::vector<int> &tour,
                        const std::vector<EndsDir> &constraints);

    /// Does \p constraints have obvious over-fixing infeasibilities.
    bool obvious_infeas(const std::vector<EndsDir> &constraints);

private:

    /// Construct a Graph::AdjList recording \p constraints.
    Graph::AdjList get_fixed_adj(const std::vector<EndsDir> &constraints);


    const Data::Instance &tsp_inst;
    Data::BestGroup &best_data;
    Graph::CoreGraph &core_graph;

    std::vector<int> fix_degrees;

    util::EdgeHash tour_edge_tracker; //!< Tracking edges added by branch tours.
    std::vector<Graph::Edge> extra_edges; //!< Extra edges for tour finding.
    int default_length; //!< Default length for sparse instances.
};

}
}

#endif
