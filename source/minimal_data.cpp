#include "datagroups.hpp"
#include "graph_io.hpp"
#include "PSEP_util.hpp"
#include "Graph.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stdexcept>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/util.h>
}

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using std::vector;
using std::unique_ptr;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace PSEP {

void Data::make_cut_test(const string &tsp_fname,
			 const string &tour_nodes_fname,
			 const string &lp_sol_fname, GraphGroup &graph_data,
			 BestGroup &best_data, vector<double> &lp_edges,
			 SupportGroup &supp_data)
{
  unique_ptr<Data::Instance> inst_p;
  Data::make_cut_test(tsp_fname, tour_nodes_fname, lp_sol_fname, graph_data,
		      best_data, lp_edges, supp_data, inst_p);
}

void Data::make_cut_test(const string &tsp_fname,
			 const string &tour_nodes_fname,
			 const string &lp_sol_fname, GraphGroup &graph_data,
			 BestGroup &best_data, vector<double> &lp_edges,
			 SupportGroup &supp_data,
			 std::unique_ptr<Data::Instance> &inst_p)
{
  int ncount = 0;
  runtime_error err("Problem in make_cut_test");

  if(inst_p)
    throw logic_error("Passed non-null Instance ptr to make_cut_test.");

  try {
    inst_p = PSEP::make_unique<Data::Instance>(tsp_fname, ncount);
  } catch(const exception &e){
    cerr << e.what() << "\n";
    throw err;
  }

  Data::Instance &inst = *inst_p;

  graph_data.m_graph.node_count = ncount;
  best_data.min_tour_value = 0.0;
  
  try {
    graph_data.island.resize(ncount);
    graph_data.edge_marks.resize(ncount, 0);
    best_data.perm.resize(ncount);
  } catch (const exception &e) {
    cerr << e.what() << " trying ncount resizes\n";
    throw err;
  }

  try {
    get_tour_nodes(ncount, best_data.best_tour_nodes, tour_nodes_fname);
    get_lp_sol(ncount, supp_data.support_elist, supp_data.support_ecap,
	       lp_sol_fname);
  } catch(const exception &e) { cerr << e.what() << "\n"; throw err; }
  
  for(int i = 0; i < ncount; ++i){//add the best tour edges
    best_data.perm[best_data.best_tour_nodes[i]] = i;
    int
      end0 = fmin(best_data.best_tour_nodes[i],
		 best_data.best_tour_nodes[(i + 1) % ncount]),
      end1 = fmax(best_data.best_tour_nodes[i],
		 best_data.best_tour_nodes[(i + 1) % ncount]);
    Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, inst.ptr()));
    best_data.min_tour_value += e.len;

    Graph &graph = graph_data.m_graph;

    try {
      graph.edges.push_back(e);
      graph.edge_lookup[IntPair(end0, end1)] = graph.edges.size() - 1;
      graph.edge_count += 1;
      best_data.best_tour_edges.push_back(1);
      graph_data.delta.push_back(0);
    } catch (const exception &e) {
      cerr << e.what() << " pushing back new tour edge. \n";
      throw err;
    }
  }

  try {
    lp_edges.resize(graph_data.m_graph.edge_count, 0);
  } catch (const exception &e) {
    cerr << e.what() << " resizing lp edges\n";
    throw err;
  }

  //add the lp edges, growing graph_data if necessary
  for(int i = 0; i < supp_data.support_ecap.size(); ++i){
    vector<int> &sup_elist = supp_data.support_elist;
    vector<double> &sup_ecap = supp_data.support_ecap;
    int end0 = sup_elist[2 * i], end1 = sup_elist[(2 * i) + 1];
    //get_lp_sol guarantees end0 < end1
    Graph &graph = graph_data.m_graph;

    IntPairMap::iterator
      edge_it = graph.edge_lookup.find(IntPair(end0, end1));

    if(edge_it == graph.edge_lookup.end()){
      Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, inst.ptr()));

      try {
	graph.edges.push_back(e);
	lp_edges.push_back(sup_ecap[i]);
	graph.edge_lookup[IntPair(end0, end1)] = graph.edges.size() - 1;
	graph.edge_count += 1;
	best_data.best_tour_edges.push_back(0);
	graph_data.delta.push_back(0);
      } catch (const exception &e) {
	cerr << e.what() << " pushing back new lp edge\n";
	throw err;
      }
    } else {
      lp_edges[edge_it->second] = sup_ecap[i];
    }
  }

  try {
    for(int i = 0; i < lp_edges.size(); ++i)
      if(lp_edges[i] >= Epsilon::Zero)
	supp_data.support_indices.push_back(i);
  } catch(const exception &e) {
    cerr << e.what() << " pushing back sup inds\n";
    throw err;
  }

  if(GraphUtils::build_s_graph(ncount, graph_data.m_graph.edges,
				   supp_data.support_indices,
			       lp_edges, &supp_data.G_s))
    throw err;
}

}
