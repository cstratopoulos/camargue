/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Data group structures.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_DATAGROUP_H
#define CMR_DATAGROUP_H

#include "util.hpp"
#include "err_util.hpp"
#include "graph.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/util.h>
}


namespace CMR {

/** Namespace for storing data groups. */
namespace Data {

/** Class for storing raw TSP instance data.
 * This class is a handle to the Concorde CCdatagroup structure, modeling
 * unique_ptr-style sole ownership over a datagroup.
 */
class Instance {
public:
    Instance() noexcept; //!< Default construct an empty instance.

    /// Construct an Instance from a TSPLIB file with a random seed.
    Instance(const std::string &fname, int seed);

    /// Construct a geometric random Instance.
    Instance(int seed, int ncount, int gridsize);

    /// Construct a sparse Instance with fixed edge set.
    Instance(const std::string &probname, int seed, int ncount,
             std::vector<int> &elist, std::vector<int> &elen);

    Instance(Instance &&I) noexcept;
    Instance &operator=(Instance &&I) noexcept;

    Instance(const Instance &I) = delete;
    Instance &operator=(const Instance &I) = delete;

    ~Instance();


    /// Access the raw pointer to the data, for use by Concorde routines.
    CCdatagroup* ptr() const { return const_cast<CCdatagroup*>(&dat); }

    /// Edge length between two nodes in an Instance.
    double edgelen(int end0, int end1) const
        { return CCutil_dat_edgelen(end0, end1, ptr()); }

    /// A function object for the edge length.
    const std::function<double(int, int)> edgelen_func() const
        { return [this](int e0, int e1){ return edgelen(e0, e1); }; }

    int node_count() const { return nodecount; } //!< Number of nodes.
    int seed() const { return random_seed; } //!< Random seed used.

    /// Get the length of the tour in \p tour_nodes.
    double tour_length(const std::vector<int> &tour_nodes) const;

    /// The TSPLIB instance name or the random problem dimensions.
    const std::string &problem_name() const { return pname; }

private:
    CCdatagroup dat; //!< The Concorde data structure being managed.

    int nodecount;
    int random_seed;

    std::string pname;
};

}

namespace Graph {

/// Edge generation protocol to use.
enum class EdgePlan {
    Linkern, //!< 10 LK tours, with quadnearest for tiny instances.
    Delaunay, //!< Delaunay triangulation.
};

class CoreGraph {
public:
    CoreGraph() = default;

    /// Construct a CoreGraph specified by a TSP instance and edge plan.
    CoreGraph(const Data::Instance &inst, Graph::EdgePlan edge_plan);

    /// Construct a CoreGraph as above with EdgePlan::Linkern edges.
    CoreGraph(const Data::Instance &inst)
        : CoreGraph(inst, EdgePlan::Linkern) {}

    /// Construct a CoreGraph from a length func and a c-array node-node elist.
    CoreGraph(int ncount, int ecount, const int *elist,
              const std::function<double(int, int)> edgelen);

    /// Construct a CoreGraph containing the nodes of a TSP tour.
    CoreGraph(const std::vector<int> &tour_nodes,
              const std::function<double(int, int)> edgelen);

    int node_count() const { return nodecount; }
    int edge_count() const { return edges.size(); }

    /// Find the index of an edge by its endpoints, returning -1 if not found.
    int find_edge_ind(int end0, int end1) const;

    Edge get_edge(int index) const { return edges[index]; }

    std::vector<Edge> &get_edges() { return edges; }
    const std::vector<Edge> &get_edges() const { return edges; }

    const AdjList &get_adj() const { return adj_list; }

    void add_edge(int end0, int end1, int len);
    void add_edge(const Edge &e);

    void remove_edges();

    template<typename numtype>
    void tour_edge_vec(const std::vector<int> &tour_nodes,
                       std::vector<numtype> &tour_edges,
                       double &tour_val) const;

private:
    std::vector<Edge> edges;
    AdjList adj_list;
    int nodecount;
};

}

namespace Data {

/** Stores information about the current best tour.
 * This structure manages the incumbent best tour for the TSP instance, as
 * well as some extra data about the tour and to facilitate computations
 * with the tour nodes and edges
 */
struct BestGroup {
    BestGroup() : min_tour_value(0) {}

    /** The LK constructor. */
    BestGroup(const Instance &inst, Graph::CoreGraph &core_graph);

    /** The tour from file constructor. */
    BestGroup(const Instance &inst, Graph::CoreGraph &core_graph,
              const std::string &tourfile);

  std::vector<int> best_tour_edges; /**< Binary vector indicating edges used
				     * in tour */
  std::vector<int> best_tour_nodes; /**< The sequence of nodes of the tour */
  std::vector<int> perm; /**< Defined by the relation
			  * `perm[best_tour_nodes[i]] = i` */

  double min_tour_value;
};

struct SupportGroup {
    SupportGroup() = default;

    SupportGroup(const std::vector<Graph::Edge> &edges,
                 const std::vector<double> &lp_x,
                 std::vector<int> &island,
                 int ncount);

    SupportGroup(SupportGroup &&SG) noexcept;

    SupportGroup &operator=(SupportGroup &&SG) noexcept;



    bool in_subtour_poly();

    std::vector<double> lp_vec;
    std::vector<int> support_indices;
    std::vector<int> support_elist;
    std::vector<double> support_ecap;

    Graph::AdjList supp_graph;

    bool connected;
    bool integral;
};

/// Load just enough data to test separation routines.
void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
                   Graph::CoreGraph &core_graph,
		   Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   Data::SupportGroup &supp_data);

/// As above, but with access to the generated Instance.
void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
                   Graph::CoreGraph &core_graph,
		   Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   Data::SupportGroup &supp_data,
		   Data::Instance &inst);

}

/////////////////////// TEMPLATE METHOD IMPLEMENTATIONS //////////////////////

namespace Graph {

template<typename numtype>
void CoreGraph::tour_edge_vec(const std::vector<int> &tour_nodes,
                              std::vector<numtype> &tour_edges,
                              double &tour_val) const
{
    using std::cerr;
    using std::endl;
    using std::runtime_error;
    using std::vector;

    tour_edges.resize(edges.size());
    std::fill(tour_edges.begin(), tour_edges.end(), 0);
    tour_val = 0.0;

    int ncount = nodecount;

    for (int i = 0; i < ncount; ++i) {
        EndPts e(tour_nodes[i], tour_nodes[(i + 1) % ncount]);
        int ind = find_edge_ind(e.end[0], e.end[1]);
        if (ind == -1) {
            cerr << e << " not in CoreGraph" << endl;
            throw runtime_error("Missing tour edge in CoreGraph");
        }
        tour_edges[ind] = 1;
        tour_val += edges[ind].len;
    }
}

}
}

#endif
