/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file                
 * @brief DATA GROUP STRUCTURE DEFINITIONS
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_DATAGROUP_H
#define PSEP_DATAGROUP_H

#include "PSEP_util.hpp"
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
    delete dat;    
  }
};

}

namespace PSEP {

/** Namespace for storing data groups.
 * Classes in this namespace store data compartmentalized to specific aspects
 * of the TSP solution process. The idea is that higher level controller
 * objects in the solution process (such as TSPSolver, PureCut, CutControl, 
 * ABC, LP::Core, etc. should be initialized with some or all of these 
 * datagroups and, if applicable, individual members of relevant groups should
 * be passed to members or methods. 
 */
namespace Data {

class Instance {
public:
  Instance() = default;
  Instance(const std::string &fname, int &ncount);
  Instance(const int seed, const int ncount, const int gridsize);

  CCdatagroup* ptr() { return handle.get(); }
  
private:
  std::unique_ptr<CCdatagroup> handle;
};

/** GraphGroup stores pure combinatorial information about the problem.
 * This structure essentially encodes a weighted graph which represents
 * the TSP instance under consideration, generally with a subset of edges
 * from the complete graph.
 */
struct GraphGroup {
  GraphGroup() = default;
  GraphGroup(const std::string &fname, std::string &probname,
	     PSEP::RandProb &randprob,
	     std::unique_ptr<CCdatagroup> &dat,
	     const bool sparse, const int quadnearest,
	     const bool dump_xy); /**< Standard constructor, @see TSPSolver. */

  
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

  /** The Lin-Kernighan tour constructor.
   * In this constructor, a tour on \p graph is constructed via a call to
   * Concorde's implementation of chained Lin-Kernighan, using the problem
   * data in \p dat. This routine may add additional edges to \p graph as 
   * needed, which also changes the size of \p delta. True values of 
   * \p save_tour and \p save_tour_edges indicate the best tour nodes and edges
   * will be saved to file. 
   * @param[in] user_seed is the random seed to be passed to the LK function,
   * to allow for reproducibility. 
   * @pre \p graph, \p delta, \p dat are all initialized by GraphGroup
   * @post `if (save_tour)` then the tour will be saved to `probname.sol`
   * @post `if (save_tour_edges)` then the tour edges will be saved to 
   * `probname_tour.x`
   * @post If the tour contains edges not originally in \p graph, these will
   * be added to \p graph, with \p delta resized accordingly. 
   * @post If `graph.node_count` is even, an extra edge will be added for
   * use in construction of an initial basis, as per Padberg-Hong (1980), 
   * with adjustments to \p graph and \p delta as above
   * @post `best_tour_edges` has length graph.edge_count and `best_tour_nodes` 
   * and `perm` have length graph.node_count
   * @post min_tour_value is the length of the tour
   */
  BestGroup(PSEP::Graph &graph, std::vector<int> &delta,
	    std::unique_ptr<CCdatagroup> &dat, const std::string &probname,
	    const int user_seed, const bool save_tour,
	    const bool save_tour_edges);
  
  /** The from-file tour constructor.
   * This constructor initializes a BestGroup using a tour specified in 
   * \p tourfile. All other operations and pre/post are the same as above. 
   * @pre \p tourfile is a cyclic permutation of the numbers 0 to 
   * `graph.node_count`
   */
  BestGroup(const std::string &tourfile,
	    PSEP::Graph &graph, std::vector<int> &delta,
	    std::unique_ptr<CCdatagroup> &dat, const std::string &probname,
	    const bool write_tour, const bool write_tour_edges);

  std::vector<int> best_tour_edges; /**< Binary vector indicating edges used
				     * in tour */
  std::vector<int> best_tour_nodes; /**< The sequence of nodes of the tour */
  std::vector<int> perm; /**< Defined by the relation 
			  * `perm[best_tour_nodes[i]] = i` */

  double min_tour_value;
};

/* This group stores objects related to the LP solver/LP relaxation */
struct LPGroup {
  LPGroup(const Graph &m_graph, PSEP::LP::Prefs &_prefs,
	  const std::vector<int> &perm);
  ~LPGroup(){PSEPlp_free(&m_lp);}

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
   * prefs - see PSEP_util.h for info
   */
  PSEPlp m_lp;  
  std::vector<double> m_lp_edges;
  std::vector<int> old_colstat;
  std::vector<int> old_rowstat;
  std::vector<int> frac_colstat;
  std::vector<int> frac_rowstat;
  PSEP::LP::Prefs prefs;
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
 * @pre \p tour_nodes_fname is as in PSEP::get_tour_nodes.
 * @pre \p lp_sol_fname is as in PSEP::get_lp_sol.
 */
void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
		   PSEP::Data::GraphGroup &graph_data,
		   PSEP::Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   PSEP::Data::SupportGroup &supp_data);

void make_cut_test(const std::string &tsp_fname,
		   const std::string &tour_nodes_fname,
		   const std::string &lp_sol_fname,
		   PSEP::Data::GraphGroup &graph_data,
		   PSEP::Data::BestGroup &best_data,
		   std::vector<double> &lp_edges,
		   PSEP::Data::SupportGroup &supp_data,
		   std::unique_ptr<PSEP::Data::Instance> &inst_p);

}
}

#endif
