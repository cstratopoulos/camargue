#include "cuts.hpp"

#include <map>
#include <utility>
#include <iostream>

#include <cmath>

using std::vector;
using std::map;
using std::pair;

using std::cout;
using std::cerr;

namespace CMR {
  
template<>
void CutQueue<HyperGraph>::push_front(const HyperGraph &H)
{
    cut_q.push_front(H);
    if (cut_q.size() > q_capacity) {
        cut_q.back().delete_refs();
        cut_q.pop_back();
    }
}

template<>
void CutQueue<HyperGraph>::push_back(const HyperGraph &H)
{
    if (cut_q.size() >= q_capacity) {
        cut_q.back().delete_refs();
        cut_q.pop_back();
    }
    cut_q.push_back(H);
}

template<>
void CutQueue<HyperGraph>::pop_front()
{
    cut_q.front().delete_refs();
    cut_q.pop_front();
}

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
                node_marks[perm[edges[j].end[1]]])
                coeff_map[j] += 1.0;
        }
    }

    rmatind.reserve(coeff_map.size());
    rmatval.reserve(coeff_map.size());

        for(pair<const int, double> &kv : coeff_map) {
            rmatind.push_back(kv.first);
            rmatval.push_back(kv.second);
        }
    
}

int CutTranslate::get_sparse_row(const HyperGraph &H, vector<int> &rmatind,
				 vector<double> &rmatval, char &sense,
				 double &rhs)
{
    int rval = 0;
    sense = 'G';

    map<int, double> coef_map;
    vector<int> body_nodes;
    int deltacount = 0;
    rmatind.clear();
    rmatval.clear();

  
    rval = H.source_setbank->extract_nodelist(*(H.set_refs[0]),
                                              body_nodes);
    CMR_CHECK_RVAL(rval, "Couldn't extract nodelist, ");
  

  
    GraphUtils::get_delta(body_nodes, edges, &deltacount, delta, edge_marks);
    rval = (deltacount == 0);
    CMR_CHECK_RVAL(rval, "Body nodes gave empty delta, cut type: "
                   << ((H.cut_type == HyperGraph::CutType::Segment) ?
                       "segment" : "blossom. "));
  

    try {
        rmatind.resize(deltacount);
        rmatval.resize(deltacount);
        for (int i = 0; i < deltacount; i++)
            coef_map[delta[i]] = 1.0;
    }
    catch (const std::bad_alloc &) {
        rval = 1; CMR_GOTO_CLEANUP("Out of memory for sparse row, ");
    }

    switch (H.cut_type) {
    case HyperGraph::CutType::Segment:    
        rhs = 2;
        break;

    case HyperGraph::CutType::Blossom:
        int num_teeth = H.set_refs.size() - 1;
        rhs = 1 - num_teeth;
    
        //skip the first ref, it is body above
        for (int i = 1; i < H.set_refs.size(); i++) {
            vector<int> edge_tooth;
            int edge_index;
            IntPairMap::iterator find_it;
      
            rval = H.source_setbank->extract_nodelist(*(H.set_refs[i]),
                                                      edge_tooth);
            CMR_CHECK_RVAL(rval, "Couldn't extract blossom tooth, ");

            if (edge_tooth.size() != 2) {
                rval = 1; CMR_GOTO_CLEANUP("Blossom tooth has "
                                           << edge_tooth.size()
                                           << "nodes! ");
            }

      
            find_it = edge_lookup.find(IntPair(fmin(edge_tooth[0],
                                                    edge_tooth[1]),
                                               fmax(edge_tooth[0],
                                                    edge_tooth[1])));
            rval = (find_it == edge_lookup.end());
            CMR_CHECK_RVAL(rval, "Couldn't find tooth in edge lookup, ");
            
            edge_index = find_it->second;
            coef_map[edge_index] = -1.0;
        }
    
        break;
    }

    {//scoped temporary variable
        int i = 0;
        for (map<int, double>::iterator it = coef_map.begin();
             it != coef_map.end(); ++it) {
            rmatind[i] = it->first;
            rmatval[i] = it->second;
            i++;
        }
    }

 CLEANUP:
    if (rval) {
        cerr << "CutTranslate<HyperGraph>::get_sparse_row failed, "
             << "row is invalid\n";
        rmatind.clear();
        rmatval.clear();
    }
    return rval;
}

int CutTranslate::get_sparse_row(const dominoparity &dp_cut,
				 const vector<int> &tour_nodes,
				 vector<int> &rmatind, vector<double> &rmatval,
				 char &sense, double &rhs)
{
    int rval = 0;
    vector<double> coeff_buff;

    for (int &i : edge_marks) i = 0;
    rmatind.clear();
    rmatval.clear();
    rhs = 0.0;
    sense = 'L';

    try { coeff_buff.resize(edges.size(), 0.0); } catch (...) {
        CMR_SET_GOTO(rval, "Couldn't allocate coefficient buffer. ");
    }

    for (const int node : dp_cut.degree_nodes)
        edge_marks[tour_nodes[node]] = 1;

    for (int i = 0; i < edges.size(); ++i) {
        Edge e = edges[i];
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
        int end0 = fmin(tour_nodes[ends.first], tour_nodes[ends.second]);
        int end1 = fmax(tour_nodes[ends.first], tour_nodes[ends.second]);

        IntPairMap::const_iterator it = edge_lookup.find(IntPair(end0, end1));
    
        if (it == edge_lookup.end()) {
            CMR_SET_GOTO(rval, "Tried to lookup invalid edge. ");
        }

        coeff_buff[it->second] -= 1.0;
    }

    try {
        for (int i = 0; i < coeff_buff.size(); ++i)
            if (coeff_buff[i] != 0.0) {
                rmatind.push_back(i);
                rmatval.push_back(coeff_buff[i]);
            }
    } catch (...) { CMR_SET_GOTO(rval,
                                 "Couldn't get sparse row from buffer. "); }

    rhs /= 2;
    rhs = floor(rhs);
  
    for (double &coeff : rmatval) {
        if (fabs(coeff >= Epsilon::Cut)) {
            coeff /= 2;
            coeff = floor(coeff);
        }
    }

 CLEANUP:
    if (rval) {
        cerr << "Problem in CutTranslate::get_sparse_row(dp_cut, ...), "
             << "cut invalidated.\n";
        rmatind.clear();
        rmatval.clear();
    }
    return rval;
}

int CutTranslate::get_sparse_row_if(bool &violated, const HyperGraph &H,
				    const vector<double> &x,
				    vector<int> &rmatind,
				    vector<double> &rmatval, char &sense,
				    double &rhs)
{
    int rval = 0;
    violated = false;
    double activity;

    rval = get_sparse_row(H, rmatind, rmatval, sense, rhs);
    CMR_CHECK_RVAL(rval, "Couldn't get sparse row, ");

    get_activity(activity, x, rmatind, rmatval);

    switch (sense) {
    case 'G':
        violated = (activity < rhs) && ((rhs - activity) >= 0.001);
        break;
    case 'L':
        violated = (activity > rhs) && ((activity - rhs) >= 0.001);
        break;
    default:
        rval = 1;
        CMR_GOTO_CLEANUP("Uncaught row sense " << sense << ", ");
    }

 CLEANUP:
    if (rval)
        cerr << "CutTranslate<HyperGraph>::get_sparse_row_if failed, "
             << "row is invalid.\n";

    if (rval || !violated) {
        rmatind.clear();
        rmatval.clear();
    }
  
    return rval;
}

int CutTranslate::is_cut_violated(bool &violated, const HyperGraph &H,
				  vector<double> &x)
{
    int rval = 0;
    vector<int> rmatind;
    vector<double> rmatval;
    char sense; double rhs;

    violated = false;

    rval = get_sparse_row_if(violated, H, x, rmatind, rmatval, sense, rhs);
    CMR_CHECK_RVAL(rval, "Couldn't test violation, ");

 CLEANUP:
    if (rval)
        cerr << "CutTranslate::is_cut_violated failed\n";
    return rval;
}

}
