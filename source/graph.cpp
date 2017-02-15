#include "graph.hpp"

#include <algorithm>
#include <fstream>		   

using std::vector;
using std::cout;
using std::cerr;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {

int Graph::connected(SupportGraph *G, int *icount,
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

void Graph::dfs(int n, SupportGraph *G, int *icount,
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

void Graph::get_delta (const std::vector<int> &nodelist,
                       std::vector<Graph::Edge> &elist,
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

void Graph::get_delta (int nsize, int *nlist, int ecount, int *elist,
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

void Graph::get_delta(const int interval_start, const int interval_end,
			   const vector<int> &tour_nodes,
			   const vector<int> &elist,
			   int &deltacount, vector<int> &delta,
			   vector<int> &node_marks)
{
  for (int i = interval_start; i <= interval_end; i++)
    node_marks[tour_nodes[i]] = 1;

  int k = 0;

  for (int i = 0; i < elist.size() / 2; i++) {
    if (node_marks[elist[2 * i]] + node_marks[elist[(2 * i) + 1]] == 1)
      delta[k++] = i;
  }

  for (int i = interval_start; i <= interval_end; i++)
    node_marks[tour_nodes[i]] = 0;  
}

int Graph::build_s_graph (int node_count,
                          const vector<Graph::Edge> &edges,
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

namespace Graph {

vector<int> delta_inds(const vector<int> &node_list, const vector<Edge> &edges,
                       int ncount)
{
    vector<int> result;
    vector<bool> node_marks(ncount, false);

    for (int n : node_list)
        node_marks[n] = true;

    for (int i = 0; i < edges.size(); ++i) {
        Edge e = edges[i];
        if (node_marks[e.end[0]] != node_marks[e.end[1]])
            result.push_back(i);
    }

    return result;
}

TourGraph::TourGraph() noexcept { CCtsp_init_lpgraph_struct(&L); }

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

TourGraph::TourGraph(TourGraph &&T) noexcept : d_tour(std::move(T.d_tour))
{
    CCtsp_free_lpgraph(&L);
    L = T.L;

    CCtsp_init_lpgraph_struct(&T.L);
    T.d_tour.clear();
}

TourGraph &TourGraph::operator=(TourGraph &&T) noexcept
{
    d_tour = std::move(T.d_tour);
    
    CCtsp_free_lpgraph(&L);
    L = T.L;

    CCtsp_init_lpgraph_struct(&T.L);
    T.d_tour.clear();

    return *this;
}


TourGraph::~TourGraph() {  CCtsp_free_lpgraph(&L); }

/**
 * This constructor is meant to be used to create an AdjList representation of
 * a CoreGraph.
 */
AdjList::AdjList(int ncount, const vector<Edge> &ref_elist) try
    : node_count(ncount), edge_count(ref_elist.size()),
      nodelist(vector<Node>(node_count))
{
    for (int i = 0; i < edge_count; ++i) {
        const Edge &e = ref_elist[i];

        nodelist[e.end[0]].neighbors.emplace_back(e.end[1], i, e.len);
        nodelist[e.end[1]].neighbors.emplace_back(e.end[0], i, e.len);
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("AdjList elist constructor failed.");
}

/**
 * Creates an AdjList representation of a support graph for an LP solution 
 * relative to some CoreGraph.
 * @param[in] ncount the number of nodes.
 * @param[in] ref_elist the Edge set of the reference CoreGraph.
 * @param[in] edge_caps capacities to be assigned to all the edges in 
 * \p ref_elist.
 * @param[in] keep_indices indices of all the edges from \p ref_elist which
 * should make it in to the AdjList.
 */
AdjList::AdjList(int  ncount,
                 const vector<Edge> &ref_elist,
                 const vector<double> &edge_caps,
                 const std::vector<int> &keep_indices) try
    : node_count(ncount), edge_count(keep_indices.size()),
      nodelist(vector<Node>(node_count))
{
    if (edge_caps.size() != ref_elist.size())
        throw logic_error("Size mismatch in support graph AdjList constructor");
    
    for (int index : keep_indices) {
        const Edge &e = ref_elist[index];
        double cap = edge_caps[index];

        nodelist[e.end[0]].neighbors.emplace_back(e.end[1], index, cap);
        nodelist[e.end[1]].neighbors.emplace_back(e.end[0], index, cap);
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("AdjList indices/ecap constructor failed.");
}

AdjList::AdjList(AdjList &&AL) noexcept
    : node_count(AL.node_count), edge_count(AL.edge_count),
      nodelist(std::move(AL.nodelist))
{}

AdjList &AdjList::operator=(Graph::AdjList &&AL) noexcept
{
    node_count = AL.node_count;
    edge_count = AL.edge_count;
    nodelist = std::move(AL.nodelist);

    return *this;
}

bool AdjList::connected(vector<int> &island, int start_node)
{
    island.clear();
    for (Node &n : nodelist)
        n.mark = 0;

    dfs(0, island);

    return island.size() == node_count;
}

void AdjList::dfs(int start_node, std::vector<int> &island)
{
    Node &n = nodelist[start_node];
    n.mark = 1;
    island.push_back(start_node);

    for (AdjObj &a : n.neighbors)
        if (nodelist[a.other_end].mark == 0)
            dfs(a.other_end, island);    
}

/**
 * The edge is added iff it is not already present.
 * @param[in] index the index, relative to some CoreGraph reference Edge set,
 * of the edge to be added.
 * @param[in] val the value which shall become the `val` field of the added
 * AdjObj.
 */
void AdjList::add_edge(int end0, int end1, int index, double val)
{
    if (find_edge(end0, end1) != nullptr)
        return;

    nodelist[end0].neighbors.emplace_back(end1, index, val);
    nodelist[end1].neighbors.emplace_back(end0, index, val);
    ++edge_count;
}

CoreGraph::CoreGraph(int ncount, int ecount, const int *elist,
                     const std::function<double(int, int)> edgelen) try
    : nodecount(ncount)
{
    edges.reserve(ecount);
    
    for (int i = 0; i < ecount; ++i) {
        int e0 = elist[2 * i];
        int e1 = elist[(2 * i) + 1];
        edges.emplace_back(e0, e1, edgelen(e0, e1));
    }

    adj_list = AdjList(ncount, edges);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("CoreGraph constructor failed.");
}

CoreGraph::CoreGraph(const vector<int> &tour_nodes,
                     const std::function<double(int, int)> edgelen) try
    : nodecount(tour_nodes.size())
{
    edges.reserve(nodecount);

    for (int i = 0; i < nodecount; ++i) {
        int e0 = tour_nodes[i];
        int e1 = tour_nodes[(i + 1) % nodecount];

        edges.emplace_back(e0, e1, edgelen(e0, e1));
    }

    adj_list = AdjList(nodecount, edges);    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("CoreGraph constructor failed.");
}

int CoreGraph::find_edge_ind(int end0, int end1) const
{
    const AdjObj *adj_ptr = adj_list.find_edge(end0, end1);
    if (adj_ptr == nullptr)
        return -1;
    return adj_ptr->edge_index;
}

void CoreGraph::get_elist(vector<int> &elist, vector<int> &elen) const
{
    elen.resize(edges.size());
    elist.resize(2 * edges.size());

    for (int i = 0; i < edges.size(); ++i) {
        const Edge &e = edges[i];
        elen[i] = e.len;
        elist[2 * i] = e.end[0];
        elist [(2 * i) + 1] = e.end[1];
    }
}

void CoreGraph::add_edge( int end0, int end1, int len )
{
    if (find_edge_ind(end0, end1) != -1)
        return;

    int new_ind = edge_count();

    edges.emplace_back(end0, end1, len);
    adj_list.add_edge(end0, end1, new_ind, len);
}

void CoreGraph::add_edge(const Edge &e)
{
    if (find_edge_ind(e.end[0], e.end[1]) != -1)
        return;

    int new_ind = edge_count();
    edges.push_back(e);
    adj_list.add_edge(e.end[0], e.end[1], new_ind, e.len);
}

}
}
