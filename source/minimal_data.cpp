#include "datagroups.hpp"
#include "graph_io.hpp"
#include "PSEP_util.hpp"
#include "Graph.hpp"

#include <vector>
#include <string>

#include <cmath>

extern "C" {
#include <concorde/INCLUDE/util.h>
}

using std::cout;
using std::cerr;
using std::vector;
using std::string;

namespace PSEP {

int Data::make_cut_test(const string &tsp_fname, const string &tour_nodes_fname,
			const string &lp_sol_fname, GraphGroup &graph_data,
			BestGroup &best_data, vector<double> &lp_edges,
			SupportGroup &supp_data)
{
  int rval = 0;
  int ncount = 0;
  CCdatagroup dat;
  best_data.min_tour_value = 0.0;

  CCutil_init_datagroup(&dat);

  rval = CCutil_gettsplib(const_cast<char *>(tsp_fname.c_str()), &ncount, &dat);
  PSEP_CHECK_RVAL(rval, "CCutil_gettsplib failed. ");

  graph_data.m_graph.node_count = ncount;
  cout << "Set ncount to " << graph_data.m_graph.node_count << ", ncount "
       << ncount << "\n";

  try {
    graph_data.island.resize(ncount);
    graph_data.edge_marks.resize(ncount, 0);
    best_data.perm.resize(ncount);
  } catch (...) {
    PSEP_SET_GOTO(rval, "Out of memory for ncount resizes. ");
  }

  rval = get_tour_nodes(ncount, best_data.best_tour_nodes, tour_nodes_fname);
  if(rval) goto CLEANUP;

  rval = get_lp_sol(ncount, supp_data.support_elist, supp_data.support_ecap,
		    lp_sol_fname);
  if(rval) goto CLEANUP;

  
  for(int i = 0; i < ncount; ++i){//add the best tour edges
    best_data.perm[best_data.best_tour_nodes[i]] = i;
    int
      end0 = fmin(best_data.best_tour_nodes[i],
		 best_data.best_tour_nodes[(i + 1) % ncount]),
      end1 = fmax(best_data.best_tour_nodes[i],
		 best_data.best_tour_nodes[(i + 1) % ncount]);
    Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, &dat));
    best_data.min_tour_value += e.len;

    Graph &graph = graph_data.m_graph;

    try {
      graph.edges.push_back(e);
      graph.edge_lookup[IntPair(end0, end1)] = graph.edges.size() - 1;
      graph.edge_count += 1;
      best_data.best_tour_edges.push_back(1);
      graph_data.delta.push_back(0);
    } catch (...) {
      PSEP_SET_GOTO(rval, "Problem pushing back new tour edge. ");
    }
  }

  try { lp_edges.resize(graph_data.m_graph.edge_count, 0); } catch (...) {
    PSEP_SET_GOTO(rval, "Couldn't make lp edges size match. ");
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
      Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, &dat));

      try {
	graph.edges.push_back(e);
	lp_edges.push_back(sup_ecap[i]);
	graph.edge_lookup[IntPair(end0, end1)] = graph.edges.size() - 1;
	graph.edge_count += 1;
	best_data.best_tour_edges.push_back(0);
	graph_data.delta.push_back(0);
      } catch (...) {
	PSEP_SET_GOTO(rval, "Problem pushing back new lp edge. ");
      }
    } else {
      lp_edges[edge_it->second] = sup_ecap[i];
    }
  }

  {
    vector<int> e_lengths = {lp_edges.size(), best_data.best_tour_edges.size(),
			     graph_data.delta.size()};
    vector<int> ref_ecount(3, graph_data.m_graph.edge_count);
    if(e_lengths != ref_ecount){
      PSEP_SET_GOTO(rval, "Edge length vector sizes do not match. ");
    }
  }

  try {
    for(int i = 0; i < lp_edges.size(); ++i)
      if(lp_edges[i] >= LP::EPSILON)
	supp_data.support_indices.push_back(i);
  } catch(...) {
    PSEP_SET_GOTO(rval, "Couldn't push back sup inds. ");
  }

  cout << "Calling build_s_graph with ncount " << ncount
       << ", ecount " << graph_data.m_graph.edge_count << "\n";
  rval = GraphUtils::build_s_graph(ncount, graph_data.m_graph.edge_count,
				   graph_data.m_graph.edges,
				   supp_data.support_indices,
				   lp_edges, &supp_data.G_s);
  if(rval) goto CLEANUP;
    

 CLEANUP:
  if(rval)
    cerr << "Problem in Data::make_cut_test.\n";
  CCutil_freedatagroup(&dat);
  return rval;
}

}
