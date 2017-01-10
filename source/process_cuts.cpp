#include "process_cuts.hpp"
#include "err_util.hpp"

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

namespace CMR {
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
    const vector<Edge> &edges = core_graph.get_edges();

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

    for (int &i : edge_marks) i = 0;
    rmatind.clear();
    rmatval.clear();
    rhs = 0.0;
    sense = 'L';

    const vector<Edge> &edges = core_graph.get_edges();

    try {
        coeff_buff.resize(edges.size(), 0.0);
    } CMR_CATCH_PRINT_THROW("allocating coeff buffer", err);

    for (const int node : dp_cut.degree_nodes)
        edge_marks[tour_nodes[node]] = 1;

    for (int i = 0; i < edges.size(); ++i) {
        const Edge &e = edges[i];
        int sum = edge_marks[e.end[0]] + edge_marks[e.end[1]];
        coeff_buff[i] += sum;
    }

    rhs += 2 * dp_cut.degree_nodes.size();

    for (const int node : dp_cut.degree_nodes)
        edge_marks[tour_nodes[node]] = 0;

    for (const SimpleTooth &T : dp_cut.used_teeth) {
        for (int i = T.body_start; i <= T.body_end; ++i)
            edge_marks[tour_nodes[i]] = 1;
        edge_marks[tour_nodes[T.root]] = -2;

        for (int i = 0; i < edges.size(); ++i) {
            Edge e = edges[i];
            int sum = edge_marks[e.end[0]] + edge_marks[e.end[1]];
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
            edge_marks[tour_nodes[i]] = 0;
        edge_marks[tour_nodes[T.root]] = 0; 
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

}
}
