#include "process_cuts.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <utility>
#include <stdexcept>

#include <cmath>

using std::vector;
using std::map;
using std::pair;

using std::cout;
using std::cerr;

using std::exception;
using std::runtime_error;
using std::logic_error;

namespace CMR {

namespace Eps = Epsilon;

namespace Sep {

void CutTranslate::get_sparse_row(const CCtsp_lpcut_in &cc_cut,
                                  const std::vector<int> &perm,
                                  vector<int> &rmatind,
                                  vector<double> &rmatval, char &sense,
                                  double &rhs)
{
    rmatind.clear();
    rmatval.clear();
    sense = cc_cut.sense;
    rhs = cc_cut.rhs;

    int ncount = perm.size();
    map<int, double> coeff_map;
    const vector<Graph::Edge> &edges = core_graph.get_edges();

    for (int i = 0; i < cc_cut.cliquecount; ++i) {
        vector<bool> node_marks(ncount, false);
        CCtsp_lpclique &clq = cc_cut.cliques[i];

        for (int j = 0; j < clq.segcount; ++j) {
            CCtsp_segment &seg = clq.nodes[j];

            for (int k = seg.lo; k <= seg.hi; ++k)
                node_marks[k] = true;
        }

        for (int j = 0; j < edges.size(); ++j) {
            if (node_marks[perm[edges[j].end[0]]] !=
                node_marks[perm[edges[j].end[1]]]) {
                if (coeff_map.count(j))
                    coeff_map[j] += 1.0;
                else
                    coeff_map[j] = 1.0;
            }
        }
    }

    rmatind.reserve(coeff_map.size());
    rmatval.reserve(coeff_map.size());

        for(pair<const int, double> &kv : coeff_map) {
            rmatind.push_back(kv.first);
            rmatval.push_back(kv.second);
        }
    
}


void CutTranslate::get_sparse_row(const dominoparity &dp_cut,
                                      const vector<int> &tour_nodes,
                                      vector<int> &rmatind,
                                      vector<double> &rmatval,
                                      char &sense, double &rhs)
{
    runtime_error err("Problem in CutTranslate::get_sparse_row dp cut.");
    vector<double> coeff_buff;

    for (int &i : node_marks) i = 0;
    rmatind.clear();
    rmatval.clear();
    rhs = 0.0;
    sense = 'L';

    const vector<Graph::Edge> &edges = core_graph.get_edges();

    try {
        coeff_buff.resize(edges.size(), 0.0);
    } CMR_CATCH_PRINT_THROW("allocating coeff buffer", err);

    for (const int node : dp_cut.degree_nodes)
        node_marks[tour_nodes[node]] = 1;

    for (int i = 0; i < edges.size(); ++i) {
        const Graph::Edge &e = edges[i];
        int sum = node_marks[e.end[0]] + node_marks[e.end[1]];
        coeff_buff[i] += sum;
    }

    rhs += 2 * dp_cut.degree_nodes.size();

    for (const int node : dp_cut.degree_nodes)
        node_marks[tour_nodes[node]] = 0;

    for (const SimpleTooth &T : dp_cut.used_teeth) {
        for (int i = T.body_start; i <= T.body_end; ++i)
            node_marks[tour_nodes[i]] = 1;
        node_marks[tour_nodes[T.root]] = -2;

        for (int i = 0; i < edges.size(); ++i) {
            Graph::Edge e = edges[i];
            int sum = node_marks[e.end[0]] + node_marks[e.end[1]];
            switch (sum) {
            case 2:
                coeff_buff[i] += 2;
                break;
            case -1:
                coeff_buff[i] += 1;
                break;
            default:
                break;
            }
        }

        rhs += (2 * (T.body_end - T.body_start + 1)) - 1;
    
        for (int i = T.body_start; i <= T.body_end; ++i)
            node_marks[tour_nodes[i]] = 0;
        node_marks[tour_nodes[T.root]] = 0; 
    }

    for (const IntPair &ends : dp_cut.nonneg_edges) {
        int e0 = tour_nodes[ends.first];
        int e1 = tour_nodes[ends.second];

        int find_ind = core_graph.find_edge_ind(e0, e1);
        if (find_ind == -1) {
            cerr << "Tried to lookup invalid edge.\n";
            throw err;
        }

        coeff_buff[find_ind] -= 1.0;
    }


    try {
        for (int i = 0; i < coeff_buff.size(); ++i)
            if (coeff_buff[i] != 0.0) {
                rmatind.push_back(i);
                rmatval.push_back(coeff_buff[i]);
            }
    } CMR_CATCH_PRINT_THROW("getting sparse row from buffer", err);

    rhs /= 2;
    rhs = floor(rhs);
  
    for (double &coeff : rmatval) {
        if (fabs(coeff >= Epsilon::Zero)) {
            coeff /= 2;
            coeff = floor(coeff);
        }
    }
}

/**
 * @param[in] B the blossom to expand.
 * @param[in] tour_edges the active tour vector.
 * @param[in] lp_vec the lp solution used to find \p B.
 * @param[in] edges the CoreGraph edges.
 * @param[in] ncount the number of nodes.
 * @param[in] handle_delta the delta_inds for `B.handle`.
 * Let \f$ \bar x \f$ be the edge vector of \p tour_edges, \f$ E^* \f$ the 
 * edges for which \p lp_vec is nonzero, and \f$ E^0, E^1 \f$ be the edges in
 * \f$ E^* \f$ for which \f$ \bar x \f$ is respectively zero or one. Let 
 * \f$ e \f$ be `B.cut_ind` and \f$ H \f$ be `B.handle`.
 * As per Letchford and Lodi (Primal Separation Algorithms), this function 
 * returns \f[ e \cup \delta(H)\cap E^1 \f] if \f$ e \in E^0 \f$, and 
 * \f[ \delta(H)\cap E^1 \setminus e \f] if \f$ e \in E^1 \f$.
 */
vector<int> teeth_inds(const ex_blossom &B, const vector<int> &tour_edges,
                       const vector<double> &lp_vec,
                       const vector<Graph::Edge> &edges, int ncount,
                       const vector<int> &handle_delta)
{
    vector<int> result = handle_delta;

    int cut_ind = B.cut_edge;

    // base teeth are delta(handle) intersect edges equal to one in tour.
    result.erase(std::remove_if(result.begin(), result.end(),
                                [&lp_vec, &tour_edges](int ind)
                                {
                                    return (lp_vec[ind] < Eps::Zero ||
                                            tour_edges[ind] != 1);
                                }),
                 result.end());

    auto it = std::find(result.begin(), result.end(), cut_ind);
    int tour_entry = tour_edges[cut_ind];

    if (tour_entry == 0) { // then we add the cut edge
        if (it == result.end())
            result.push_back(tour_entry);
    } else if (tour_entry == 1) { // then we remove it
        if (it != result.end())
            result.erase(it);
    } else
        throw logic_error("Non-binary tour entry");

    return result;
}


/// As above, but without precomputed handle_
vector<int> teeth_inds(const ex_blossom &B, const vector<int> &tour_edges,
                       const vector<double> &lp_vec,
                       const vector<Graph::Edge> &edges, int ncount)
{
    const vector<int> &handle = B.handle;
    
    vector<int> handle_delta = Graph::delta_inds(handle, edges, ncount);

    return teeth_inds(B, tour_edges, lp_vec, edges, ncount, handle_delta);
}

/**
 * Same arguments as teeth_inds. Returns true if `B.handle` is too small
 * or too big, if there are an even number of teeth returned by teeth_inds,
 * or if the teeth intersect.
 */
bool bad_blossom(const ex_blossom &B, const vector<int> &tour_edges,
                 const vector<double> &lp_vec,
                 const vector<Graph::Edge> &edges, int ncount)
{
    const vector<int> &handle = B.handle;
    
    if (handle.size() < 3 || handle.size() > ncount - 3)
        return true;

    vector<int> teeth = teeth_inds(B, tour_edges, lp_vec, edges, ncount);
    
    if ((teeth.size() % 2) == 0)
        return true;

    vector<int> node_marks(ncount, 0);

    //now check for intersecting teeth.
    for (int ind : teeth) {
        const Graph::Edge &e = edges[ind];
        for (int n : e.end) {
            ++node_marks[n];
            if (node_marks[n] > 1)
                return true;
        }
    }

    return false;    
}

}
}
