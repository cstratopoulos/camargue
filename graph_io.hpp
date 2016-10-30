/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *
 *                UTILITIES FOR READING AND WRITING TO FILES
 *
 * This file contains some simple functions for file input/output related to
 * tours and edge sets.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_GRAPH_IO_HPP
#define PSEP_GRAPH_IO_HPP

#include "Graph.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace PSEP {
/*
 * Writes the tour specified by tour_nodes to the file named tour_nodes_fname
 * PRE:
 * tour_nodes should be nonempty, a cyclic permutation of the numbers 0 to 
 *            tour_nodes.size()
 * tour_nodes_fname should be a nonempty string
 * POST:
 * tour_nodes_fname will be the name of a file with tour_nodes.size() at the
 *            top, followed by the entries of tour_nodes
 * returns 0 if successful, 1 if error
 */
int write_tour_nodes(const std::vector<int> &tour_nodes,
		     const std::string &tour_nodes_fname);

/*
 * Writes the edges in edges specified by tour_edges to the file named
 * tour_edges_fname, assuming the graph has node_count nodes
 * PRE:
 * tour_edges and edges are nonempty, of the same size
 * tour_edges is binary, and the collection of edges[i] for which 
 * tour_edges[i] == 1 gives a connected cyclical graph
 * edges is a list of edges in a graph with node_count nodes
 * tour_edges_fname is a nonempty string
 * POST:
 * tour_edges_fname will be the name of a file with
 * node_count edges.size()
 * at the top, followed by
 * edges[i].end[0] edges[i].end[1] 1.0
 * for each i such that tour_edges[i] = 1.0
 * returns 0 if successful, 1 if error
 */
int write_tour_edges(const std::vector<int> &tour_edges,
		     const std::vector<PSEP::Edge> &edges,
		     const int node_count,
		     const std::string &tour_edges_fname);

/*
 * Writes the tour specified by lp_elist and lp_ecap to the file named
 * lp_edges_fname
 * PRE:
 * lp_elist and lp_ecap nonempty,
 * lp_elist.size() is twice lp_ecap.size()
 * lp_ecap[i] is the weight on the edge lp_elist[2i], lp_elist[2i + 1]
 * lp_edges_fname is a nonempty string
 * POST:
 * lp_edges_fname will be the name of a file with 
 * node_count lp_ecap.size() 
 * at the top, followed by 
 * support_elist[2i] support_elist[2i + 1] support_ecap[2i]
 * for all i from 0 to lp_ecap.size()
 * returns 0 if successful, 1 if error
 */
int write_lp_edges(const std::vector<int> &lp_elist,
		     const std::vector<double> &lp_ecap,
		     const std::string &lp_edges_fname);

/*
 * Stores the tour specified in tour_nodes_fname to the vector tour_nodes
 * PRE: 
 * tour_nodes_fname names an existant, nonempty file with a node count at 
 * the top, followed by a cyclic permutation of the numbers from 0 to that 
 * node_count
 * POST:
 * tour_nodes has length specified by the first line of tour_nodes_fname,
 * and stores its subsequent entries in the same order
 * returns 0 if successful, 1 if error
 */
int get_tour_nodes(std::vector<int> &tour_nodes,
		   const std::string &tour_nodes_fname);
		     
}


#endif
