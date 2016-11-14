#include "graph_io.hpp"
#include "PSEP_util.hpp"

#include <algorithm>
#include <sstream>

#include <cstdio>
#include <cmath>


namespace PSEP {

int write_tour_nodes(const std::vector<int> &tour_nodes,
		     const std::string &tour_nodes_fname)
{
  int rval = 0;
  std::ofstream tour_out;

  if(tour_nodes.empty())
    PSEP_SET_GOTO(rval, "Tried to write empty tour to file. ");

  if(tour_nodes_fname.empty())
    PSEP_SET_GOTO(rval, "Tried to specfy empty filename. ");

  try {
    tour_out.open(tour_nodes_fname);
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't open outfile stream. ");
  }
  
  try {
    tour_out << tour_nodes.size() << "\n";

    int i;
    for(i = 0; i < tour_nodes.size(); ++i){
      tour_out << tour_nodes[i] << " ";
      if(i % 10 == 9)
	tour_out << "\n";
    }

    if(i % 10)
      tour_out << "\n";    
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't write tour to file. ");
  }

 CLEANUP:
  tour_out.close();
  if(rval){
    std::cerr << "PSEP::write_tour_nodes failed.\n";
    std::remove(tour_nodes_fname.c_str());
  }
  return rval;
}

int write_tour_edges(const std::vector<int> &tour_edges,
		     const std::vector<PSEP::Edge> &edges,
		     const int node_count,
		     const std::string &tour_edges_fname)
{
  int rval = 0;
  int ecount = 0;
  std::ofstream tour_out;

  for(int i : tour_edges)
    ecount += i;

  if(tour_edges.empty() || edges.empty())
    PSEP_SET_GOTO(rval, "Tried to write empty tour to file. ");

  if(tour_edges_fname.empty())
    PSEP_SET_GOTO(rval, "Tried to specify empty filename. ");

  if(tour_edges.size() != edges.size())
    PSEP_SET_GOTO(rval, "Sizes of edges and tour_edges incompatible. ");

  try {
    tour_out.open(tour_edges_fname);
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't open outfile stream. ");
  }

  try {
    tour_out << node_count << " " << ecount << "\n";

    for(int i = 0; i < tour_edges.size(); ++i){
      if(tour_edges[i] == 1)
	tour_out << edges[i].end[0] << " " << edges[i].end[1] << " 1.0\n";
    }
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't write tour to file. ");
  }

 CLEANUP:
  tour_out.close();
  if(rval){
    std::cerr << "PSEP::write_tour_edges failed.\n";
    std::remove(tour_edges_fname.c_str());
  }
  return rval;
}

int write_lp_edges(const std::vector<int> &lp_elist,
		   const std::vector<double> &lp_ecap,
		   const int node_count,
		   const std::string &lp_edges_fname)
{
  int rval = 0;
  std::ofstream lp_out;

  if(node_count <= 0)
    PSEP_SET_GOTO(rval, "Passed bad value of node count. ");

  if(lp_elist.empty() || lp_ecap.empty())
    PSEP_SET_GOTO(rval, "Passed empty lp solution. ");

  if(lp_edges_fname.empty())
    PSEP_SET_GOTO(rval, "Tried to write to empty filename. ");

  if(lp_elist.size() != 2 * lp_ecap.size())
    PSEP_SET_GOTO(rval, "Incompatible elist and ecap. ");

  try {
    lp_out.open(lp_edges_fname);
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't open outfile stream. ");
  }

  try {
    lp_out << node_count << " " << lp_ecap.size() << "\n";

    for(int i = 0; i < lp_ecap.size(); ++i){
      lp_out << lp_elist[2 * i] << " " << lp_elist[(2 * i) + 1] << " "
	     << lp_ecap[i] << "\n";
    }
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't write LP solution to file. ");
  }

 CLEANUP:
  lp_out.close();
  if(rval){
    std::cerr << "PSEP::write_lp_edges failed\n";
    std::remove(lp_edges_fname.c_str());
  }
  return rval;
}

int write_xy_coords(const double *x, const double *y, const int node_count,
		    const std::string &xy_coords_fname)
{
  int rval = 0;
  std::ofstream xy_out;

  if(node_count <= 0)
    PSEP_SET_GOTO(rval, "Passed bad value of ncount. ");

  if(!x || !y)
    PSEP_SET_GOTO(rval, "Passed null x or y coords. ");

  if(xy_coords_fname.empty())
    PSEP_SET_GOTO(rval, "Tried to specify empty filename. ");

  try {
    xy_out.open(xy_coords_fname);
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't open outfile stream. ");
  }

  try {
    xy_out << node_count << "\n";

    for(int i = 0; i < node_count; ++i)
      xy_out << x[i] << " " << y[i] << "\n";
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't write xy coords. ");
  }

 CLEANUP:
  xy_out.close();
  if(rval){
    std::cerr << "PSEP::write_xy_coords failed\n";
    std::remove(xy_coords_fname.c_str());
  }
  return rval;
}

int get_tour_nodes(const int node_count, std::vector<int> &tour_nodes,
		   const std::string &tour_nodes_fname)
{
  int rval = 0;
  std::ifstream tour_in;
  std::string current_line;
  std::vector<int> uniq_vec;
  int file_ncount;

  if(node_count == 0)
    PSEP_SET_GOTO(rval, "Specified zero nodes. ");

  if(tour_nodes_fname.empty())
    PSEP_SET_GOTO(rval, "Specified empty filename. ");

  try { tour_in.open(tour_nodes_fname); } catch (std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't open infile stream. ");
  }

  if(!tour_in.good())
    PSEP_SET_GOTO(rval, "Bad infile. ");

  try { tour_in >> file_ncount; } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't read in file nodecount. ");
  }

  if(file_ncount != node_count)
    PSEP_SET_GOTO(rval, "File node count does not match specified count. ");

  try {
    tour_nodes.reserve(node_count);
    uniq_vec.reserve(node_count);
  } catch(std::bad_alloc &e) {
    PSEP_SET_GOTO(rval, "Couldn't reserve tour nodes or uniq vec. ");
  }

  try {
    while(std::getline(tour_in, current_line)){
      std::stringstream line_stream(current_line);
      int current_node;

      while(line_stream >> current_node)
	tour_nodes.push_back(current_node);
    }
  } catch(std::ios_base::failure &e) {
    PSEP_SET_GOTO(rval, "Couldn't read in tour nodes. ");
  }

  if(tour_nodes.size() != node_count)
    PSEP_SET_GOTO(rval, "File contains wrong number of nodes. ");

  uniq_vec = tour_nodes;

  std::sort(uniq_vec.begin(), uniq_vec.end());

  if (std::unique(uniq_vec.begin(), uniq_vec.end()) != uniq_vec.end())
    PSEP_SET_GOTO(rval, "Tour input file contains duplicate nodes. ");

  if(uniq_vec.front() != 0 || uniq_vec.back() != node_count - 1)
    PSEP_SET_GOTO(rval, "Tour input file contains nodes out of range. ");
  

 CLEANUP:
  tour_in.close();
  if(rval){
    std::cerr << "PSEP::get_tour_nodes failed\n";
    tour_nodes.clear();
  }
  return rval;
}

int get_lp_sol(const int node_count, std::vector<int> &support_elist,
	       std::vector<double> &support_ecap,
	       const std::string &lp_sol_fname)
{
  int rval = 0;
  std::ifstream lp_in;
  int file_ncount, file_ecount;

  if(node_count == 0) PSEP_SET_GOTO(rval, "Specified zero nodes. ");

  if(lp_sol_fname.empty()) PSEP_SET_GOTO(rval, "Gave empty filename. ");

  try { lp_in.open(lp_sol_fname); } catch (std::ios_base::failure &) {
    PSEP_SET_GOTO(rval, "Couldn't open infile stream. ");
  }

  if(!lp_in.good()) PSEP_SET_GOTO(rval, "Infile is not good. ");

  try { lp_in >> file_ncount >> file_ecount; }
  catch (std::ios_base::failure &) {
    PSEP_SET_GOTO(rval, "Couldn't read in file ncount or ecount. ");
  }

  if(file_ncount != node_count) PSEP_SET_GOTO(rval, "File has wrong ncount. ");

  try {
    support_ecap.resize(file_ecount);
    support_elist.resize(2 * file_ecount);
  } catch (std::bad_alloc &) {
    PSEP_SET_GOTO(rval, "Out of memory for sup vecs. ");
  }

  try {
    for(int i = 0; i < file_ecount; ++i){
      int end0, end1;
      lp_in >> end0 >> end1 >> support_ecap[i];
      support_elist[2 * i] = fmin(end0, end1);
      support_elist[(2 * i) + 1] = fmax(end0, end1);
    }
  } catch (std::ios_base::failure &) {
    PSEP_SET_GOTO(rval, "Problem reading file lines. ");
  }

 CLEANUP:
  if(rval)
    std::cerr << "get_lp_sol failed.\n";
  lp_in.close();
  return rval;
}


}
