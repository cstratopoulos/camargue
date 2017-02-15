/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file                
 * @brief Data group structures. 
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_DATAGROUP_H
#define CMR_DATAGROUP_H

#include "util.hpp"
#include "graph.hpp"

#include <functional>
#include <memory>
#include <string>
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

    /// The TSPLIB instance name or the random problem dimensions.
    const std::string &problem_name() const { return pname; }
  
private:
    CCdatagroup dat; //!< The Concorde data structure being managed.

    int nodecount;
    int random_seed;
    
    std::string pname;
};

/** Pure combinatorial information about the problem.
 * This structure uses a CoreGraph to encode a weighted undirected graph 
 * containing the subset of complete graph edges under consideration in the
 * current CoreLP.
 */
struct GraphGroup {
    GraphGroup() = default;

    /// Generate edges from an Instance and create the associated Coregraph.
    GraphGroup(const Instance &inst);

    Graph::CoreGraph core_graph; //!< The Edge list and AdjList.
    
    std::vector<int> island; //!< Stores components from a dfs of m_graph.
    std::vector<int> delta; //!< Stores edges in delta of some node set.
    std::vector<int> node_marks; //!< Marks nodes for adjacency computations.
};

/** Stores information about the current best tour.
 * This structure manages the incumbent best tour for the TSP instance, as
 * well as some extra data about the tour and to facilitate computations
 * with the tour nodes and edges
 */
struct BestGroup {
    BestGroup() : min_tour_value(0) {}

    /** The LK constructor. */
    BestGroup(const Instance &inst, GraphGroup &graph_data);

    /** The tour from file constructor. */
    BestGroup(const Instance &inst, GraphGroup &graph_data,
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

/** Load just enough Data to test cut separation routines.
 * This pseudo-constructor function is designed for testing separation routines
 * using the tour in \p tour_nodes_fname and the lp solution in \p lp_sol_fname
 * as a tour/lp solution on the tsp instance in \p tsp_fname. It will 
 * populate \p graph_data with precisely the edges in \p lp_sol_fname, together
 * with the edges joining adjacent nodes in \p tour_nodes_fname. Then, 
 * \p best_data will be initialized with the tour in \p tour_nodes_fname. 
 * The vector \p lp_edges will then have size equal to
 * `graph_data.m_graph.edge_count`, with binary entries indicating the lp
 * solution from \p lp_sol_fname. Finally, \p lp_edges and \p graph_data
 * will be used to populate \p supp_data. 
 * @pre \p tsp_fname is the name of a TSPLIB instance.
 * @pre \p tour_nodes_fname is as in CMR::get_tour_nodes.
 * @pre \p lp_sol_fname is as in CMR::get_lp_sol.
 */
void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
		   Data::GraphGroup &graph_data,
		   Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   Data::SupportGroup &supp_data);

void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
		   Data::GraphGroup &graph_data,
		   Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   Data::SupportGroup &supp_data,
		   Data::Instance &inst_p);

}
}

#endif
