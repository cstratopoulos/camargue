#include "simpleDP.hpp"
#include "DPgraph.hpp"
#include "PSEP_util.hpp"

#include <memory>

#include <cmath>

using std::vector;
using std::cout;
using std::cerr;

namespace PSEP {


int Cut<dominoparity>::separate()
{
  int rval = 0;
  std::unique_ptr<DPCutGraph> witness;
  
  rval = candidates.get_light_teeth();
  if(rval) goto CLEANUP;

  candidates.weak_elim();

  witness = PSEP::make_unique<DPCutGraph>(candidates);
  rval = witness->simple_DP_sep(dp_q);

 CLEANUP:
  if(rval == 1)
    cerr << "Cut<dominoparity>::separate failed.\n";
  return rval;
}

int Cut<dominoparity>::parse_cut(const dominoparity &dp_cut,
				 Data::GraphGroup &g_dat, SupportGraph &G_s,
				 vector<double> &rmatval, double &rhs)
{
  int rval = 0;
  
  if(rmatval.size() != g_dat.m_graph.edge_count)
    PSEP_SET_GOTO(rval, "rmatval is not the right size! ");

  for(const SimpleTooth *T : dp_cut.used_teeth)
    parse_tooth(T, g_dat, G_s, rmatval, rhs);

  parse_handle(dp_cut.degree_nodes, g_dat, rmatval, rhs);

  for(const IntPair &ends : dp_cut.nonneg_edges)
    parse_nonneg_edges(ends, g_dat, rmatval);

  for(double &coeff : rmatval)
    coeff /= 2;

  rhs /= 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Cut<dominoparity>::parse_cuts failed.\n";
  return rval;
}

void Cut<dominoparity>::parse_tooth(const SimpleTooth *T,
				    Data::GraphGroup &g_dat,
				    SupportGraph &G_s,
				    vector<double> &rmatval, double &rhs)
{
  vector<int> &node_marks = g_dat.edge_marks;
  vector<Edge> &edges = g_dat.m_graph.edges;

  //adding 2|S| - 1.
  rhs += (2 * (T->body_end - T->body_start + 1)) - 1;

  for(int i = T->body_start; i <= T->body_end; ++i)
    node_marks[best_tour_nodes[i]] = 1;

  //adding 2x(E(S))
  for(int i = 0; i < edges.size(); ++i)
    if(node_marks[edges[i].end[0]] + node_marks[edges[i].end[1]] == 2)
      rmatval[i] += 2;

  SNode &current_node = G_s.nodelist[best_tour_nodes[T->root]];
  for(int i = 0; i < current_node.s_degree; ++i)
    if(node_marks[current_node.adj_objs[i].other_end] == 1)
      rmatval[current_node.adj_objs[i].edge_index] += 1;

  for(int i = T->body_start; i <= T->body_end; ++i)
    node_marks[best_tour_nodes[i]] = 0;
}

void Cut<dominoparity>::parse_handle(const vector<int> &handle_nodes,
				     Data::GraphGroup &g_dat,
				     vector<double> &rmatval, double &rhs)
{
  vector<int> &node_marks = g_dat.edge_marks;
  vector<Edge> &edges = g_dat.m_graph.edges;
  
  rhs += 2 * handle_nodes.size();

  for(int i : handle_nodes)
    node_marks[best_tour_nodes[i]] = 1;

  for(int i = 0; i < edges.size(); ++i)
    rmatval[i] += node_marks[edges[i].end[0]] + node_marks[edges[i].end[1]];
  
  for(int i : handle_nodes)
    node_marks[best_tour_nodes[i]] = 0;
}

void Cut<dominoparity>::parse_nonneg_edges(const IntPair &edge_ends,
					   Data::GraphGroup &g_dat,
					   vector<double> &rmatval)
{
  int end0 = fmin(best_tour_nodes[edge_ends.first],
		  best_tour_nodes[edge_ends.second]);
  int end1 = fmax(best_tour_nodes[edge_ends.first],
		  best_tour_nodes[edge_ends.second]);
  IntPairMap::iterator it = g_dat.m_graph.edge_lookup.find(IntPair(end0,
								   end1));
  rmatval[(*it).second] -= 1;
}

int Cut<dominoparity>::cutcall()
{
  int rval = separate();

 CLEANUP:
  if(rval == 1)
    cerr << "Cut<dominoparity>::cutcall failed.\n";
  return rval;
}

}
