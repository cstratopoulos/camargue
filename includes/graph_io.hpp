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

/** Writes a tour to file.
 * Writes the tour specified by \p tour_nodes to the file named 
 * \p tour_nodes_fname
 * @pre \p tour_nodes is nonempty, containing a cyclic permutation of the 
 * numbers from 0 to \p tour_nodes.size()
 * @pre \p tour_nodes_fname is a nonempty filename
 * @post \p tour_nodes_fname will be the name of a file with 
 * tour_nodes.size() at the top, followed by the entries of tour_nodes
 * \return 0 if successful, 1 if error
 */
void write_tour_nodes(const std::vector<int> &tour_nodes,
		     const std::string &tour_nodes_fname);

/** Writes tour edges to file.
 * Writes the edges in edges specified by \p tour_edges to the file named
 * \p tour_edges_fname, assuming the graph has \p node_count nodes.
 * @param[in] edges a vector of Edge structs indicating edge endpoints
 * @pre \p tour_edges and \p edges are nonempty, of the same size
 * @pre \p tour_edges is binary, and the collection of \p edges[i] for which 
 * \p tour_edges`[i] == 1` gives a connected cyclical graph
 * @pre \p edges is a list of edges in a graph with \p node_count nodes
 * @pre \p tour_edges_fname is a nonempty string
 * @post \p tour_edges_fname will be the name of a file with \p node_count
 * and  `edges.size()` on the top line, followed by \p edges`[i].end[0]`
 * \p edges`[i].end[1]` `1.0` for each i such that \p tour_edges`[i] == 1`
 * \return 0 if successful, 1 if error
 */
void write_tour_edges(const std::vector<int> &tour_edges,
		     const std::vector<PSEP::Edge> &edges,
		     const int node_count,
		     const std::string &tour_edges_fname);

/** Writes LP edges to file.
 * Writes the tour specified by \p lp_elist and \p lp_ecap to the file named
 * \p lp_edges_fname
 * @param[in] lp_elist the edges in an LP solution in node node format
 * @param[in] lp_ecap the value assigned to each edge in the LP solution
 * @pre \p lp_elist and \p lp_ecap nonempty
 * @pre `lp_elist.size()` is twice `lp_ecap.size()`
 * @pre `lp_ecap[i]` is the weight on the edge `lp_elist[2i]` , 
 *  `lp_elist[2i + 1]`
 * @pre \p lp_edges_fname is a nonempty string
 * @post \p lp_edges_fname will be the name of a file with
 *  \p node_count `lp_ecap.size()` on the first line, followed by
 * `support_elist[2i] support_elist[2i + 1] support_ecap[2i]`
 * for all `i` from 0 to  `lp_ecap.size()`
 * \return 0 if successful, 1 if error
 */
void write_lp_edges(const std::vector<int> &lp_elist,
		   const std::vector<double> &lp_ecap,
		   const int node_count,
		   const std::string &lp_edges_fname);

/** Dumps the xy-coordinates for nodes in a graph.
 * Populates \p xy_coords_fname with the xy-coords specified by \p x, \p y
 * for a graph with \p ncount edges
 * @pre \p x and \p y are non-null double arrays of length \p ncount
 * @post \p xy_coords_fname has \p ncount on its first line, followed by
 * \p x`[i]` \p y`[i]` on all the following lines. 
 */
void write_xy_coords(const double *x, const double *y, const int ncount,
		    const std::string &xy_coords_fname);

/** Loads a tour from file.
 * Stores the tour specified in \p tour_nodes_fname to the vector
 * \p tour_nodes for a TSP instance on \p node_count cities
 * @pre \p tour_nodes_fname names an existant file whose first line is
 * \p node_count and whose following entries are a cyclic permutation of the 
 * numbers 0, ..., \p node_count
 * @post \p tour_nodes has length \p node_count and its entries are the nodes
 * from \p tour_nodes_fname in the same order
 * \returns 0 if successful, 1 if error
 */
void get_tour_nodes(const int node_count, std::vector<int> &tour_nodes,
		   const std::string &tour_nodes_fname);


/** Loads an lp solution from file.
 * Reads the lp solution specified in \p lp_sol_fname, storing the ends of
 * the edges in \p support_elist and the corresponding lp weight in 
 * \p support_ecap, for a graph with \p node_count edges. 
 * @pre \p lp_sol_fname names an existant file whose first line is `n m` where
 * `n` is \p node_count and `m` is the number of nonzero edges in the solution.
 * @post `support_ecap.size() == m` and `support_elist.size() == 2 * m`.
 *
 * @post For all positive `i`, if the `i`th line of \p lp_sol_fname is
 *  `u v w` then
 * @post `support_elist[2 * i] == u`
 * @post `suppor_elist[(2 * i) + 1] == v`
 * @post `support_ecap[i] == w`
 */
void get_lp_sol(const int node_count, std::vector<int> &support_elist,
	       std::vector<double> &support_ecap,
	       const std::string &lp_sol_fname);
		     
}


#endif
