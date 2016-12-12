/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file                
 * @brief DATA GROUP STRUCTURE DEFINITIONS
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_DATAGROUP_H
#define CMR_DATAGROUP_H

#include "util.hpp"
#include "lp.hpp"
#include "Graph.hpp"

extern "C" {
#include <concorde/INCLUDE/util.h>
}

#include <vector>
#include <memory>
#include <string>

#include <cmath>

namespace std {

template<>
struct default_delete<CCdatagroup> {
  void operator()(CCdatagroup *dat) const {
    if(dat) CCutil_freedatagroup(dat);
    CC_IFFREE(dat, CCdatagroup);
  }
};

}

namespace CMR {

/** Namespace for storing data groups.
 * Classes in this namespace store data compartmentalized to specific aspects
 * of the TSP solution process. The idea is that higher level controller
 * objects in the solution process (such as TSPSolver, PureCut, CutControl, 
 * ABC, LP::Core, etc. should be initialized with some or all of these 
 * datagroups and, if applicable, individual members of relevant groups should
 * be passed to members or methods. 
 */
namespace Data {

/** Class for storing raw TSP instance data.
 * This class is a handle to the Concorde CCdatagroup structure, modeling 
 * unique_ptr-style sole ownership over a datagroup. Move and copy operations
 * are essentially those which exist on the underlying unique_ptr.
 */
class Instance {
public:
    Instance() noexcept; /**< Default construct an instance with null handle. */

    /** Construct an instance from a TSPLIB file.
     * If \p fname is the path to a TSPLIB file from the executable directory,
     * constructs an Instance specifying the data in \p fname.
     */
    Instance(const std::string &fname, const int seed);

    /** Construct a geometric random TSP instance.
     * The instance will have \p ncount nodes distributed uniform randomly over
     * the \p gridsize by \p gridsize grid, with random seed \p seed.
     */
    Instance(const int seed, const int ncount, const int gridsize);
  
    Instance(const Instance &I) = delete;
    Instance(Instance &&I) noexcept; /**< Move constructor. */

  
    Instance &operator=(const Instance &I) = delete;
    Instance &operator=(Instance &&I) noexcept; /**< Move assign. */
  
    /** Access the raw pointer to the data, for use by Concorde routines. */
    CCdatagroup* ptr() const { return const_cast<CCdatagroup*>(handle.get()); }

    int node_count() const { return nodecount; }
    int seed() const { return random_seed; }
    const std::string &problem_name() const { return pname; }
  
private:
    std::unique_ptr<CCdatagroup> handle;

    int nodecount;
    int random_seed;
    
    std::string pname;
};

/** GraphGroup stores pure combinatorial information about the problem.
 * This structure essentially encodes a weighted graph which represents
 * the TSP instance under consideration, generally with a subset of edges
 * from the complete graph.
 */
struct GraphGroup {
    GraphGroup() = default;
    GraphGroup(const Instance &inst);
  
    Graph m_graph; /**< A Graph object describing the TSP instance. */
    
    std::vector<int> island; /**< Stores components from a dfs of m_graph */
    std::vector<int> delta; /**< Stores edges in delta of some node set */
    std::vector<int> edge_marks; /**< Marks nodes for adjacency computations */
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

/* This group stores objects related to the LP solver/LP relaxation */
struct LPGroup {
  LPGroup(const Graph &m_graph, CMR::LP::Prefs &_prefs,
	  const std::vector<int> &perm);
  ~LPGroup(){CMRlp_free(&m_lp);}

  /*
   * m_lp - the LP environment/problem object for use with the routines 
   *    in lp.h
   * m_lp_edges - vector of length graph.edge_count:
   *    the most recently computed LP solution with entries corresponding
   *    to weights assigned in the solution
   * The colstat vectors have entries equal to the symbolic constants
   * CPX_AT_LOWER, CPX_BASIC, or CPX_AT_UPPER. 
   * The 'old' ones store the basis associated with the current best tour
   * The 'frac' ones store the basis associated with the last LP solution
   * prefs - see util.h for info
   */
  CMRlp m_lp;  
  std::vector<double> m_lp_edges;
  std::vector<int> old_colstat;
  std::vector<int> old_rowstat;
  std::vector<int> frac_colstat;
  std::vector<int> frac_rowstat;
  CMR::LP::Prefs prefs;
};

/* 
 * SupportGroup is the structure responsible for managing a support graph
 * and the information about the associated LP solution
 */
struct SupportGroup  {
  /*
   * G_s - a graph whose edges are the edges from GraphGroup::m_graph
   *     for which the corresponding entry of m_lp_edges is nonnegative
   * support_indices - a list of the nonnegative edge indices
   * support_elist - the edges in support_indices, in node node format
   * support_ecap - the value assigned to the edges in support_elist in the
   *    current LP solution
   */
  SupportGraph G_s;
  std::vector<int> support_indices;
  std::vector<int> support_elist;
  std::vector<double> support_ecap;
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
		   CMR::Data::GraphGroup &graph_data,
		   CMR::Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   CMR::Data::SupportGroup &supp_data);

void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
		   CMR::Data::GraphGroup &graph_data,
		   CMR::Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   CMR::Data::SupportGroup &supp_data,
		   CMR::Data::Instance &inst_p);

}
}

#endif
