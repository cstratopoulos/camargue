/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Working with branch tours.
 */ /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_BRANCH_TOUR_H
#define CMR_BRANCH_TOUR_H

#include "branch_node.hpp"
#include "datagroups.hpp"
#include "graph.hpp"
#include "core_lp.hpp"
#include "edgehash.hpp"

#include <utility>
#include <vector>


namespace CMR {
namespace ABC {

/// Compute branch tours for estimation and instatement, managing their edges.
class BranchTourFind {
public:
    BranchTourFind(const Data::Instance &inst, const Data::BestGroup &bestdata,
                   const Graph::CoreGraph &coregraph, LP::CoreLP &corelp);

    /// Compute a tour to be instated at \p B, modifying the core LP if needed.
    void instate_branch_tour(const BranchNode &B, bool &found_tour,
                             std::vector<int> &tour);

    /// Compute a branch tour estimate, returning feasibility and tour length.
    void estimate_tour(const std::vector<EndsDir> &constraints,
                       bool &feas, double &tour_val);

    /// Is the \p tour compliant with \p constraints.
    bool tour_compliant(const std::vector<int> &tour,
                        const std::vector<EndsDir> &constraints);

    /// Does \p constraints have obvious over-fixing infeasibilities.
    bool obvious_infeas(const std::vector<EndsDir> &constraints);

    /// Call chained Lin-Kernighan to compute a branch tour.
    void compute_tour(const std::vector<EndsDir> &edge_stats,
                      bool &found_tour, bool &feas,
                      std::vector<int> &tour, double &tour_val,
                      bool for_use);

    /// Compute a list of common constraints for splitting on \p parent.
    std::vector<EndsDir> common_constraints(const BranchNode &parent,
                                            const EndPts &branch_edge);

    /// Compute a list of branch constraints for \p N and all its ancestors.
    std::vector<EndsDir> branch_constraints(const BranchNode &N);

    int verbose = 0;

private:
    /// Construct a Graph::AdjList recording \p constraints.
    Graph::AdjList get_fixed_adj(const std::vector<EndsDir> &constraints);


    const Data::Instance &tsp_inst;
    const Data::BestGroup &best_data;
    const Graph::CoreGraph &core_graph;
    LP::CoreLP &core_lp;

    std::vector<int> fix_degrees; //!< Tracking degrees for fixed up edges.

    std::vector<Graph::Edge> extra_edges; //!< Extra edges for tour finding.

    /// Default length assigned to edges not explictly in sparse instances.
    int default_length;

    /// Large length for sparse instances.
    /// This cost will be assigned to all 0-fixed edges, and `-large_length`
    /// will be assigned to all 1-fixed edges.
    int large_length;
};

}
}

#endif
