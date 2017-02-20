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
namespace Graph {

void get_elist(const vector<Edge> &edges, vector<int> &elist,
               vector<int> &elen)
{
    elen.clear();
    elist.clear();
    elen.resize(edges.size());
    elist.resize(2 * edges.size());

    for (int i = 0; i < edges.size(); ++i) {
        const Edge &e = edges[i];
        elen[i] = e.len;
        elist[2 * i] = e.end[0];
        elist [(2 * i) + 1] = e.end[1];
    }
}

void get_delta(const int interval_start, const int interval_end,
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

vector<int> delta_inds(const vector<int> &node_list, const vector<int> &elist,
                       int ncount)
{
    vector<int> result;
    vector<bool> node_marks(ncount, false);

    for (int n : node_list)
        node_marks[n] = true;

    if (elist.size() % 2 != 0)
        throw logic_error("Odd number of nodes in elist for delta_inds");

    int ecount = elist.size() / 2;

    for (int i = 0; i < ecount; ++i)
        if (node_marks[elist[2 * i]] != node_marks[elist[(2 * i) + 1]])
            result.push_back(i);

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

AdjList &AdjList::operator=(AdjList &&AL) noexcept
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

}
}
