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
 * The following structures are provided in the PSEP::Data:: namespace
 *
 * GraphGroup is the data group for holding pure combinatorial information
 * about the problem.
 * GraphGroup(string &fname, RandProb &randprob, unique_ptr<CCdatagroup> &dat)
 * See tsp_solver.h for info on initializing with fname/randprob/dat
 *
 * Members
 *                m_graph: A Graph object describing the TSP instance, see
 *                    Graph.h for information
 *                island: a vector of ints where each int represents a node
 *                    in a connected component found during a depth first 
 *                    search. Re-used repeatedly throughout execution
 *                delta: vector of ints representing a list of edges
 *                    constituting the graph theoretic notion of delta of a set
 *                    of nodes, i.e., edges with one end in the set and one
 *                    end out of it. The set is indicated by
 *                edge_marks: a vector of size equal to the number of nodes
 *                    in the graph, with an entry of one to indicate membership
 *                    in a set being used to define delta.
 *
 * BestGroup stores information about the best TSP tour that we have for the 
 * current problem.
 * BestGroup(Graph &graph, unique_ptr<CCdatagroup> &dat)
 * The constructor initializes BestGroup by computing a Lin-Kernighan tour
 * and using it to populate the other member objects
 *                graph: the Graph object for which we want to find a tour
 *                dat: A CCdatagroup pointer having been initialized by the
 *                    GraphGroup constructor with information about the problem
 *                    to be used by Concorde's Lin-Kernighan solver
 *
 * Members
 *                best_tour_edges: a vector of length equal to the number of
 *                    edges in the graph, effectively a binary vector with
 *                    a 1 indicating that the edge appears in the current best
 *                    known tour, 0 meaning it is not used.
 *                best_tour_nodes: a vector of length equal to number of nodes
 *                    in the graph, with each entry being a unique integer from
 *                    0, ..., nodecount - 1, indicating an order of traversing
 *                    the nodes that gives the current best tour
 *                perm: a permutation vector. It is defined by the mapping
 *                    perm[best_tour_nodes[i]] = i, hence giving a relabelling
 *                    of the vertices in the best tour such that the sequence
 *                    0, 1, .... , nodecount - 1 corresponds to the best known
 *                    tour. 
 *                
 *                min_tour_value: the length of the current best known tour
 *
 * LPGroup stores information about the linear programming relaxation of the 
 * problem we are considering. 
 * LPGroup(Graph &m_graph, LP::Prefs &_prefs, vector<int> &_perm)
 * Initialized with m_graph from GraphGroup, prefs from the TSPSolver 
 * constructor, and perm from BestGroup.
 *
 * Members:
 *                m_lp: The linear programming problem structure defined in
 *                    lp.h, for use with CPLEX routines
 *                m_lp_edges: a vector of length equal to the number of edges
 *                    in the graph, with [0, 1] entries corresponding to the
 *                    weight an edge is given in an LP solution. In theory this
 *                    should always be a vertex of the TSP polytope adjacent
 *                    to the current best tour, or the best tour itself. 
 *                old_colstat, old_rowstat/frac_colstat, frac_rowstat:
 *                    these are integer vectors of length
 *                    equal to the number of columns (colstat) / rows (rowstat)
 *                    in m_lp. Their entries are one of the symbolic constants
 *                    CPX_AT_LOWER, CPX_BASIC, CPX_AT_UPPER, representing the
 *                    simplex algorithm basis associated with a given solution.
 *                    The old_ vectors store the basis associated with the best
 *                    tour, and the frac_ vectors store the basis associated
 *                    with the lp solution. 
 *                    These allow us to pivot back and forth from fractional
 *                    solution to best tour, as required by separation routines
 *                    and for the control flow of the primal cutting plane 
 *                    solvers. 
 *                prefs: preferences for LP solution methods, see PSEP_util.h
 *
 * SupportGroup is used for storing information related to the support graph
 * of an LP solution, a weighted graph with edges given by the edges of the TSP
 * problem instance with nonzero value in the LP solution, and edge capacities
 * given by the weight in the LP solution. Used for various separation routines
 * At every new LP solve, the routine GraphUtils::build_s_graph is used to 
 * generate the support graph associated with the current LP solution, see 
 * Graph.h for more info
 *
 * Members:
 *                G_s: a SupportGraph structure; see Graph.h
 *                support_indices: an integer vector storing the indices of 
 *                    edges that assume a nonzero value in the current LP 
 *                    solution
 *                 support_elist: the edges in support_indices, in node-node
 *                    format. That is, the ith edge has ends
 *                     support_elist[2 * i], suppport_elist[2 * i + 1]
 *                support_ecap: support_ecap[i] indicates the value of edge
 *                    i in the associated LP solution. 
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_DATAGROUP_H
#define PSEP_DATAGROUP_H

#include<vector>
#include<memory>

#include<math.h>

#include "PSEP_util.h"
#include "Graph.h"
#include "lp.h"

namespace PSEP {
  namespace Data {
    struct GraphGroup {
      GraphGroup(const std::string &fname, PSEP::RandProb &randprob,
		 std::unique_ptr<CCdatagroup> &dat);
  
      Graph m_graph;
      std::vector<int> island;
      std::vector<int> delta;
      std::vector<int> edge_marks;
    };

    struct BestGroup {
      BestGroup(const Graph &graph, std::unique_ptr<CCdatagroup> &dat);
      std::vector<int> best_tour_edges;
      std::vector<int> best_tour_nodes;
      std::vector<int> perm;

      double min_tour_value;
    };

    struct LPGroup {
      LPGroup(const Graph &m_graph, PSEP::LP::Prefs &_prefs,
	      const std::vector<int> &perm);
      ~LPGroup(){PSEPlp_free(&m_lp);}
      PSEPlp m_lp;  
      std::vector<double> m_lp_edges;
      std::vector<int> old_colstat;
      std::vector<int> old_rowstat;
      std::vector<int> frac_colstat;
      std::vector<int> frac_rowstat;
      PSEP::LP::Prefs prefs;
    };

    struct SupportGroup {
      SupportGraph G_s;
      std::vector<int> support_indices;
      std::vector<int> support_elist;
      std::vector<double> support_ecap;
    };
  }
}

#endif
