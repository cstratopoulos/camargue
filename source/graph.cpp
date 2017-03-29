#include "graph.hpp"

#include <algorithm>
#include <fstream>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

using std::runtime_error;
using std::logic_error;
using std::exception;

namespace CMR {
namespace Graph {

/**
 * @param[in] edges the Graph edges to represent.
 * @param[out] the node-node elist representation.
 * @param[out] the edge weights.
 * See the implementation for the way edges correspond to entries in \p elist
 * and \p elen.
 */
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
    if (edge_caps.size() != ref_elist.size()) {
        cerr << "Edge caps size " << edge_caps.size() << " vs ref elist size "
             << ref_elist.size() << endl;
        throw logic_error("Mismatch in support graph AdjList constructor");
    }

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
