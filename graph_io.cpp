#include "graph_io.hpp"

#include <cstdio>

#define PSEP_SET_GOTO(rval, message) {		\
    rval = 1;					\
    std::cerr << message;			\
    goto CLEANUP; }


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
    std::cerr << "write_tour_nodes failed.\n";
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
    std::cerr << "write_tour_edges failed.\n";
    std::remove(tour_edges_fname.c_str());
  }
  return rval;
}


}

int main()
{
  
  
  return 0;
}
