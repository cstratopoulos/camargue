#include "io_util.hpp"
#include "util.hpp"

#include <iomanip>
#include <algorithm>
#include <sstream>
#include <stdexcept>

#include <cstdio>
#include <cmath>

using std::runtime_error;
using std::logic_error;
using std::exception;

using std::cerr;
using std::endl;


namespace CMR {
namespace util {

/** 
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
		     const std::string &tour_nodes_fname)
{
  std::ofstream tour_out;
  std::runtime_error err("Problem in write_tour_nodes.");

  if (tour_nodes.empty())
    throw logic_error("Tried to write empty tour to file. ");

  if (tour_nodes_fname.empty())
    throw logic_error("Tried to specfy empty filename. ");

  try {
    tour_out.open(tour_nodes_fname);
  } catch (const exception &e) {
    cerr << e.what() << "opening tour_out \n";
    throw err;
  }
  
  try {
    tour_out << tour_nodes.size() << "\n";

    int i;
    for (i = 0; i < tour_nodes.size(); ++i) {
      tour_out << tour_nodes[i] << " ";
      if (i % 10 == 9)
	tour_out << "\n";
    }

    if (i % 10)
      tour_out << "\n";    
  } catch (const exception &e) {
    cerr << e.what() << " writing tour_out\n";
    throw err;
  }
}

/** 
 * @param[in] tour_edges the edges of the tour to write to file.
 * @param[in] edges a vector of Edge structs indicating edge endpoints
 * @param[in] node_count the number of nodes in the problem, hence indicating
 * maximum index of an endpoint in \p edges or \p tour_edges.
 * @param[in] tour_edges_fname the file name to write to.
 * @pre \p tour_edges and \p edges are nonempty, of the same size
 * @pre \p tour_edges is binary, and the collection of \p edges[i] for which 
 * \p tour_edges`[i] == 1` gives a connected cyclical graph
 * @pre \p edges is a list of edges in a graph with \p node_count nodes
 * @pre \p tour_edges_fname is a nonempty string
 * @post \p tour_edges_fname will be the name of a file with \p node_count
 * and  `edges.size()` on the top line, followed by \p edges`[i].end[0]`
 * \p edges`[i].end[1]` `1.0` for each i such that \p tour_edges`[i] == 1`
 */
void write_tour_edges(const std::vector<int> &tour_edges,
		     const std::vector<Graph::Edge> &edges,
		     const int node_count,
		     const std::string &tour_edges_fname)
{
  int ecount = 0;
  std::ofstream tour_out;
  runtime_error err("Problem in write_tour_edges");

  for (int i : tour_edges)
    ecount += i;

  if (tour_edges.empty() || edges.empty())
    throw logic_error("Tried to write empty tour to file. ");

  if (tour_edges_fname.empty())
    throw logic_error("Tried to specify empty filename. ");

  if (tour_edges.size() != edges.size())
    throw logic_error("Sizes of edges and tour_edges incompatible. ");

  try {
    tour_out.open(tour_edges_fname);
  } catch (const exception &e) {
    cerr << e.what() << " opening tour out\n";
    throw err;
  }

  try {
    tour_out << node_count << " " << ecount << "\n";

    for (int i = 0; i < tour_edges.size(); ++i) {
      if (tour_edges[i] == 1)
	tour_out << edges[i].end[0] << " " << edges[i].end[1] << " 1.0\n";
    }
  } catch (const exception &e) {
    cerr << e.what() << " writing tour out\n";
    throw err;
  }
}


/**
 * Writes the specified lp solution with edges and capacities to file.
 * @param[in] lp_elist the edges in an LP solution in node node format
 * @param[in] lp_ecap the value assigned to each edge in the LP solution
 * @param[in] node_count the number of nodes in the instance.
 * @param[in] lp_edges_fname the file to write to.
 * @pre \p lp_elist and \p lp_ecap nonempty
 * @pre `lp_elist.size()` is twice `lp_ecap.size()`
 * @pre `lp_ecap[i]` is the weight on the edge `lp_elist[2i]` , 
 *  `lp_elist[2i + 1]`
 * @pre \p lp_edges_fname is a nonempty string
 * @post \p lp_edges_fname will be the name of a file with
 *  \p node_count `lp_ecap.size()` on the first line, followed by
 * `support_elist[2i] support_elist[2i + 1] support_ecap[2i]`
 * for all `i` from 0 to  `lp_ecap.size()`
 */
void write_lp_edges(const std::vector<int> &lp_elist,
		   const std::vector<double> &lp_ecap,
		   const int node_count,
		   const std::string &lp_edges_fname)
{
  std::ofstream lp_out;
  runtime_error err("Problem in write_lp_edges");

  if (node_count <= 0)
    throw logic_error("Passed bad value of node count. ");

  if (lp_elist.empty() || lp_ecap.empty())
    throw logic_error("Passed empty lp solution. ");

  if (lp_edges_fname.empty())
    throw logic_error("Tried to write to empty filename. ");

  if (lp_elist.size() != 2 * lp_ecap.size())
    throw logic_error("Incompatible elist and ecap. ");

  try {
    lp_out.open(lp_edges_fname);
  } catch (const exception &e) {
    cerr << e.what() << "opening lp_out\n";
    throw err;
  }

  try {
    lp_out << node_count << " " << lp_ecap.size() << "\n";

    for (int i = 0; i < lp_ecap.size(); ++i) {
      lp_out << lp_elist[2 * i] << " " << lp_elist[(2 * i) + 1] << " "
	     << std::fixed << std::setprecision(6) << lp_ecap[i] << "\n";
      lp_out.unsetf(std::ios_base::fixed);
    }
  } catch (const exception &e) {
    cerr << e.what() << "writing lp out\n";
    throw err;
  }
}


/**
 * @param[in] x a c-style array of x coordinates
 * @param[in] y a c-style array of y coordinates
 * @param[in] node_count the number of nodes in the instance, hence the lengths
 * of \p x and \p y.
 * @param[in] xy_coords_fname the filename to write to.
 * @post \p xy_coords_fname has \p ncount on its first line, followed by
 * \p x`[i]` \p y`[i]` on all the following lines. 
 */
void write_xy_coords(const double *x, const double *y, const int node_count,
		    const std::string &xy_coords_fname)
{
  std::ofstream xy_out;
  runtime_error err("Problem in write_xy_cords");

  if (node_count <= 0)
    throw logic_error( "Passed bad value of ncount. ");

  if (!x || !y) throw logic_error("Passed null x or y coords. ");

  if (xy_coords_fname.empty())
    throw logic_error("Tried to specify empty filename. ");

  try {
    xy_out.open(xy_coords_fname);
  } catch (const exception &e) {
    cerr << e.what() << " opening xy coords\n";
    throw err;
  }

  try {
    xy_out << node_count << "\n";

    for (int i = 0; i < node_count; ++i)
      xy_out << x[i] << " " << y[i] << "\n";
  } catch (const exception &e) {
    cerr << e.what() << "wrting xy coords\n";
    throw err;
  }
}


/**
 * Reads a tour from file, storing it in a vector of nodes. Some basic checks
 * are performed to ensure that the tour is on the right number of nodes, and
 * visits no city more than once.
 * @param[in] node_count the number of nodes in the instance.
 * @param[out] tour_nodes the vector used to store the nodes.
 * @param[in] tour_nodes_fname the filename to read from.
 * @pre \p tour_nodes_fname names an existant file whose first line is
 * \p node_count and whose following entries are a cyclic permutation of the 
 * numbers 0, ..., \p node_count
 * @post \p tour_nodes has length \p node_count and its entries are the nodes
 * from \p tour_nodes_fname in the same order
 */
void get_tour_nodes(const int node_count, std::vector<int> &tour_nodes,
		   const std::string &tour_nodes_fname)
{
    std::ifstream tour_in;
    std::string current_line;
    int file_ncount;
    runtime_error err("Problem in get_tour_nodes");

    if (node_count == 0) throw logic_error( "Specified zero nodes. ");

    if (tour_nodes_fname.empty())
        throw logic_error("Specified empty filename. ");

    try { tour_in.open(tour_nodes_fname); } catch (const exception &e) {
        cerr << e.what() << " trying to open tour_in ifstream\n";
        throw err;
    }

    if (!tour_in.good())
        throw logic_error("tour_in is not good. ");

    try { tour_in >> file_ncount; } catch (const exception &e) {
        cerr << e.what() << " trying to read tour_in nodecount\n";
        throw err;
    }

    if (file_ncount != node_count)
        throw logic_error("tour_in node count doesn't match specified count. ");

    try {
        tour_nodes.resize(node_count);
    } catch (std::bad_alloc &e) {
        cerr << e.what() << " allocating tour vectors\n";
        throw err;
    }

    try {
        int i = 0;
        while (std::getline(tour_in, current_line)) {
            std::stringstream line_stream(current_line);
            int current_node;

            while (line_stream >> current_node)
                tour_nodes[i++] = current_node;
                // tour_nodes.push_back(current_node);
        }
    } catch (const exception &e) {
        cerr << e.what() << " reading tour nodes\n"; throw err;
    }

    if (tour_nodes.size() != node_count) {
        cerr << "Tour nodes size: " << tour_nodes.size() << ", node count: "
             << node_count << "\n";
        throw logic_error("Tour nodes size mismatch. ");
    }
    std::vector<int> uniq_vec = tour_nodes;

    std::sort(uniq_vec.begin(), uniq_vec.end());

    if (std::unique(uniq_vec.begin(), uniq_vec.end()) != uniq_vec.end())
        throw logic_error("Tour input file contains duplicate nodes. ");

    if (uniq_vec.front() != 0 || uniq_vec.back() != node_count - 1)
        throw logic_error("Tour input file contains nodes out of range. ");
}


/**
 * Reads a specified solution from file, storing edge endpoints and 
 * capacities in vectors.
 * @param[in] node_count the dimension of the instance.
 * @param[out] support_elist a node-node list of the edges in the solution.
 * @param[out] support_ecap capacities corresponding to the edges in 
 * support_elist.
 * @param[in] lp_sol_fname the filename to read from.
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
	       const std::string &lp_sol_fname)
{
  std::ifstream lp_in;
  int file_ncount, file_ecount;
  runtime_error err("Problem in get_lp_sol");

  if (node_count == 0)
    throw logic_error("Zero nodes in get_lp_sol");

  if (lp_sol_fname.empty())
    throw logic_error("lp_sol_fname.empty()");

  try { lp_in.open(lp_sol_fname); } catch (const exception &e) {
    cerr << e.what() << " opening lp_sol infile.";
  }

  if (!lp_in.good()) throw logic_error("lp_in is not good.");

  try { lp_in >> file_ncount >> file_ecount; } catch (const exception &e) {
    cerr << e.what() << " reading ncount or ecount\n";
    throw err;
  }

  if (file_ncount != node_count)
    throw logic_error("lp_in has wrong ncount.");
  
  try {
    support_ecap.resize(file_ecount);
    support_elist.resize(2 * file_ecount);
  } catch (const exception &e) {
    cerr << e.what() << " resizing sup vecs\n";
    throw err;
  }

  try {
    for (int i = 0; i < file_ecount; ++i) {
      int end0, end1;
      lp_in >> end0 >> end1 >> support_ecap[i];
      support_elist[2 * i] = fmin(end0, end1);
      support_elist[(2 * i) + 1] = fmax(end0, end1);
    }
  } catch (const exception &e) {
    cerr << e.what() << " reading support elist/ecap\n";
    throw err;
  }
}

}
}
