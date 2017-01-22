/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ /**
 * @file
 * @brief Utilities for reading and writing to files.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef CMR_GRAPH_IO_HPP
#define CMR_GRAPH_IO_HPP

#include "graph.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

namespace CMR {
namespace util {

/// Writes a tour to file.
void write_tour_nodes(const std::vector<int> &tour_nodes,
		     const std::string &tour_nodes_fname);

/// Writes tour edges to file.
void write_tour_edges(const std::vector<int> &tour_edges,
		     const std::vector<Graph::Edge> &edges,
		     const int node_count,
		     const std::string &tour_edges_fname);

/// Writes LP edges to file.
void write_lp_edges(const std::vector<int> &lp_elist,
		   const std::vector<double> &lp_ecap,
		   const int node_count,
		   const std::string &lp_edges_fname);

/// Dumps the xy-coordinates for nodes in a graph.
void write_xy_coords(const double *x, const double *y, const int ncount,
		    const std::string &xy_coords_fname);

/// Load a tour from file.
void get_tour_nodes(const int node_count, std::vector<int> &tour_nodes,
		   const std::string &tour_nodes_fname);


/// Loads an lp solution from file.
void get_lp_sol(const int node_count, std::vector<int> &support_elist,
	       std::vector<double> &support_ecap,
	       const std::string &lp_sol_fname);

}		     
}


#endif
