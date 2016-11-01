/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                
 *                      DATA GROUP STRUCTURE DEFINITIONS
 *
 * This header file includes several different structures that are intended
 * to group data within several categories as used by modules in this program.
 * Roughly, the idea is that higher level classes (e.g., TSPSolver, PureCut,
 * ABC, CutControl LPCore) should be initialized with some or all of these
 * datagroups, and (if applicable) they should initialize their member classes
 * with individual members of the relevant data groups.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef PSEP_DATAGROUP_H
#define PSEP_DATAGROUP_H

#include "PSEP_util.hpp"
#include "lp.hpp"
#include "Graph.hpp"

#include <vector>
#include <memory>
#include <string>

#include <cmath>

namespace PSEP {
namespace Data {  
  /* GraphGroup stores pure combinatorial information about the problem */
  struct GraphGroup {
    /* constructor parameters are exactly as in tsp_solver.hpp */
    GraphGroup(const std::string &fname, std::string &probname,
	       PSEP::RandProb &randprob,
	       std::unique_ptr<CCdatagroup> &dat,
	       const bool sparse, const int quadnearest,
	       const bool dump_xy);
    /*
     * m_graph: A Graph object describing the TSP instance, see Graph.h for 
     *    information
     * island: a vector of ints where each int represents a node
     *     in a connected component found during a depth first 
     *     search. Re-used repeatedly throughout program execution
     * delta: vector of ints representing a list of edges
     *     constituting the graph theoretic notion of delta of a set
     *     of nodes, i.e., edges with one end in the set and one
     *     end out of it. The set is indicated by....
     * edge_marks: a vector of size equal to the number of nodes
     *     in the graph, with an entry of one to indicate membership
     *     in a set being used to define delta.
     */
    Graph m_graph;
    std::vector<int> island;
    std::vector<int> delta;
    std::vector<int> edge_marks;
  };

  /* Stores information about the current best tour */
  struct BestGroup {
    /*
     * Takes graph and dat initialized by GraphGroup to construct a starting
     * tour via Lin-Kernighan
     */
    BestGroup(PSEP::Graph &graph, std::vector<int> &delta,
	      std::unique_ptr<CCdatagroup> &dat, const std::string &probname,
	      const int user_seed, const bool write_tour,
	      const bool write_tour_edges);
    /*
     * Takes graph and dat initialized by GraphGroup to construct a starting
     * tour from the file specified by tourfile
     */
    BestGroup(const std::string &tourfile,
	      PSEP::Graph &graph, std::vector<int> &delta,
	      std::unique_ptr<CCdatagroup> &dat, const std::string &probname,
	      const bool write_tour, const bool write_tour_edges);

    /*
     * best_tour_edges - a binary vector of length graph.edge_count indicating
     *    which edges are used in the current best tour
     * best_tour_nodes - an integer vector of length graph.node_count giving
     *    a sequence of nodes which correspond to the current best tour
     * perm - a permutation vector of length graph.node_count, defined by the
     *    property perm[best_tour_nodes[i]] = i, hence giving a relabelling
     *    of the nodes so that 0, 1, 2, ... graph.node_count - 1 gives the 
     *    best known tour
     */
    std::vector<int> best_tour_edges;
    std::vector<int> best_tour_nodes;
    std::vector<int> perm;

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
}
}

#endif
