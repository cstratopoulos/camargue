#include "Graph.hpp"

#include <algorithm>
#include <fstream>		   
#include <iostream>
#include <stdexcept>

using std::vector;
using std::cout;
using std::cerr;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

Edge::Edge(int e0, int e1, int _len):
    len(_len),
    removable(false)
{
  end[0] = e0 < e1 ? e0 : e1;
  end[1] = e1 > e0 ? e1 : e0;
}

void Graph::print_edges() {
    for (int i = 0; i < edges.size(); ++i)
    {
        cout << edges[i].end[0] << ", " << edges[i].end[1]
	     << "-(" << edges[i].len << ")\n";
    }
}

TourGraph::TourGraph(const vector<int> &tour_edges,
		     const vector<Edge> &edges, const vector<int> &perm)
try
{
    vector<int> elist;
    int ncount = perm.size();
    int ecount = edges.size();

    for (const Edge &e : edges) {
        elist.push_back(perm[e.end[0]]);
        elist.push_back(perm[e.end[1]]);
    }
    
  for (int i : tour_edges)
      d_tour.push_back(i);

  CCtsp_init_lpgraph_struct(&L);
  
  if (CCtsp_build_lpgraph(&L, ncount, ecount, &elist[0], (int *) NULL))
      throw runtime_error("CCtsp_build_lpgraph failed.");
  
  if (CCtsp_build_lpadj(&L, 0, ecount))
      throw runtime_error("CCtsp_build_lpadj failed.");
  
} catch (const std::exception &e) {
    cerr << e.what() << "\n";
    throw std::runtime_error("TourGraph constructor failed.");
}

TourGraph::~TourGraph() { if (&L) CCtsp_free_lpgraph(&L); }

int GraphUtils::connected(SupportGraph *G, int *icount,
			  std::vector<int> &island,
			  int starting_node) {
  *icount = 0;
  for (int i = 0; i < G->node_count; i++)
    G->nodelist[i].mark = 0;

  dfs(starting_node, G, icount, island);

  if (*icount == G->node_count)
    return 1;
  else
    return 0;
}

void GraphUtils::dfs(int n, SupportGraph *G, int *icount,
		     std::vector<int> &island)
{
  int neighbor;
  SNode *pn;

  island[*icount] = n;
  (*icount)++;

  pn = &G->nodelist[n];
  pn->mark = 1;

  for (int i = 0; i < pn->s_degree; i++) {
    neighbor = pn->adj_objs[i].other_end;
    if (G->nodelist[neighbor].mark == 0)
      dfs(neighbor, G, icount, island);
  }
}

void GraphUtils::get_delta (const std::vector<int> &nodelist,
			    std::vector<Edge> &elist,
			    int *deltacount_p, std::vector<int> &delta,
			    std::vector<int> &marks) {
  for (int i = 0; i < nodelist.size(); i++)
    marks[nodelist[i]] = 1;

  int k = 0;
  for (int i = 0; i < elist.size(); i++)
    if (marks[elist[i].end[0]] + marks[elist[i].end[1]] == 1)
      delta[k++] = i;

  *deltacount_p = k;
  
  for (int i = 0; i < nodelist.size(); i++)
    marks[nodelist[i]] = 0;
}

void GraphUtils::get_delta (int nsize, int *nlist, int ecount, int *elist,
			 int *deltacount, int *delta, int *marks) {
  int i, k = 0;

  for (i = 0; i < nsize; i++) marks[nlist[i]] = 1;

  for (i = 0; i < ecount; i++) {
    if (marks[elist[2*i]] + marks[elist[2*i+1]] == 1) {
      delta[k++] = i;
    }
  }

  *deltacount = k;

  for (i = 0; i < nsize; i++) marks[nlist[i]] = 0;
}

void GraphUtils::get_delta(const int interval_start, const int interval_end,
			   const vector<int> &tour_nodes,
			   const vector<int> &elist,
			   int &deltacount, vector<int> &delta,
			   vector<int> &edge_marks)
{
  for (int i = interval_start; i <= interval_end; i++)
    edge_marks[tour_nodes[i]] = 1;

  int k = 0;

  for (int i = 0; i < elist.size() / 2; i++) {
    if (edge_marks[elist[2 * i]] + edge_marks[elist[(2 * i) + 1]] == 1)
      delta[k++] = i;
  }

  for (int i = interval_start; i <= interval_end; i++)
    edge_marks[tour_nodes[i]] = 0;  
}

int GraphUtils::build_s_graph (int node_count,
			       const vector<Edge> &edges,
			       const vector<int> &support_indices,
			       const vector<double> &m_lp_edges,
                               SupportGraph *G_s)
{
  int i, ind, a, b;
  SNode *n;
  s_adjobj *p;
  int edge_count = support_indices.size();

  if (G_s->nodelist) free(G_s->nodelist);
  if (G_s->adjspace) free(G_s->adjspace);
  G_s->nodelist = (SNode *) malloc(node_count * sizeof(SNode));
  G_s->adjspace = (s_adjobj *) malloc(2 * edge_count * sizeof(SNode));
  if (!G_s->nodelist || !G_s->adjspace) {
    fprintf(stderr, "Out of memory for support nodelist or adjspace\n");
    return 1;
  }

  for (i = 0; i < node_count; i++)
    G_s->nodelist[i].s_degree = 0;

  for (i = 0; i < edge_count; i++) {
    ind = support_indices[i];
    a = edges[ind].end[0]; b = edges[ind].end[1];
    G_s->nodelist[a].s_degree++;
    G_s->nodelist[b].s_degree++;
  }

  p = G_s->adjspace;
  for (i = 0; i < node_count; i++) {
    G_s->nodelist[i].adj_objs = p;
    p += G_s->nodelist[i].s_degree;
    G_s->nodelist[i].s_degree = 0;
  }

  for (i = 0; i < edge_count; i++) {
    ind = support_indices[i];
    a = edges[ind].end[0]; b = edges[ind].end[1];
    n = &G_s->nodelist[a];
    n->adj_objs[n->s_degree].other_end = b;
    n->adj_objs[n->s_degree].edge_index = ind;
    n->adj_objs[n->s_degree].lp_weight = m_lp_edges[ind];
    n->s_degree++;
    n = &G_s->nodelist[b];
    n->adj_objs[n->s_degree].other_end = a;
    n->adj_objs[n->s_degree].edge_index = ind;
    n->adj_objs[n->s_degree].lp_weight = m_lp_edges[ind];
    n->s_degree++;
  }

  G_s->node_count = node_count; G_s->edge_count = edge_count;

  return 0;
}

namespace GraphUtils {

AdjList::AdjList(int ncount, const vector<CMR::Edge> &ref_elist) try
    : node_count(ncount), edge_count(ref_elist.size()),
      nodelist(vector<Node>(node_count, Node((2 * edge_count) / node_count)))
{
    for (int i = 0; i < edge_count; ++i) {
        const CMR::Edge &e = ref_elist[i];

        nodelist[e.end[0]].neighbors.emplace_back(e.end[1], i, e.len);
        nodelist[e.end[1]].neighbors.emplace_back(e.end[0], i, e.len);
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("AdjList elist constructor failed.");
}

AdjList::AdjList(int  ncount,
                 const vector<CMR::Edge> &ref_elist,
                 const vector<double> &edge_caps,
                 const std::vector<int> &keep_indices) try
    : node_count(ncount), edge_count(keep_indices.size()),
      nodelist(vector<Node>(node_count, Node((2 * edge_count) / node_count)))
{
    for (int index : keep_indices) {
        const CMR::Edge &e = ref_elist[index];
        double cap = edge_caps[index];

        nodelist[e.end[0]].neighbors.emplace_back(e.end[1], index, cap);
        nodelist[e.end[1]].neighbors.emplace_back(e.end[0], index, cap);
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("AdjList indices/ecap constructor failed.");
}

void AdjList::add_edge(int end0, int end1, int index, double val)
{
    if (find_edge(end0, end1) != nullptr)
        return;

    nodelist[end0].neighbors.emplace_back(end1, index, val);
    nodelist[end1].neighbors.emplace_back(end0, index, val);
    ++edge_count;
}

}

}
