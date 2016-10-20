#include "tooth.hpp"

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <iostream>
#include <algorithm>

using namespace PSEP;
using namespace std;

int CandidateTeeth::dump_cut(double cut_val, int cut_start, int cut_end,
			     void *u_data)
{
  linsub_cb_data *cb_data = (linsub_cb_data *) u_data;
  seg *old_cut = cb_data->old_cut;
  vector<vector<SimpleTooth::Ptr>> &lite_t = cb_data->cb_lite_teeth;
  double slack = ((double) cut_val - 2.0)/2.0;

  // for(int i = cut_start; i <= cut_end; i++)
  //   cout << cb_data->cb_tour_nodes[i] << ", ";
  // cout << "\n";

  if(cut_start == old_cut->start && cut_end == old_cut->end + 1 &&
     slack + old_cut->cutval < .4999){
    lite_t[cut_end].emplace_back(SimpleTooth::Ptr(new SimpleTooth(cut_end, cut_start, old_cut->end, slack + old_cut->cutval)));
  }

  old_cut->start = cut_start;
  old_cut->end = cut_end;
  old_cut->cutval = slack;

  return 0;
}

CandidateTeeth::CandidateTeeth(vector<int> &_edge_marks,
			       vector<int> &_best_tour_nodes,
			       SupportGraph &_G_s,
			       vector<int> &_support_elist,
			       vector<double> &_support_ecap) :
  edge_marks(_edge_marks), best_tour_nodes(_best_tour_nodes), G_s(_G_s),
  support_elist(_support_elist), support_ecap(_support_ecap)
{
  light_teeth.resize(best_tour_nodes.size());
}

int CandidateTeeth::get_light_teeth()
{
  int rval = 0;

  vector<int> endmark(G_s.node_count, CC_LINSUB_BOTH_END);
  seg lin_seg(G_s.node_count -1, G_s.node_count - 1, 0);

  linsub_cb_data cb_data = {light_teeth, best_tour_nodes, &lin_seg};

  clear_collection();

  //  print_vec(best_tour_nodes);
  
  cout << "Getting root adjacent light teeth via linsub...." << endl;
  double st = zeit();
  rval = CCcut_linsub_allcuts(G_s.node_count,
			      G_s.edge_count,
			      &best_tour_nodes[0], &endmark[0],
			      &support_elist[0], &support_ecap[0], 2.999,
			      &cb_data, dump_cut);
  if(rval) goto CLEANUP;
  st = zeit() - st;

  cout << st << "s total\n";


 CLEANUP:
  if(rval == 1){
    cerr << "CandidateTeeth::get_light_teeth failed\n";
    clear_collection();
  }
  return rval;
}

int CandidateTeeth::get_adjacent_teeth(const int root)
{
  int rval = 0;
  int ncount = G_s.node_count;
  int body_start = (root + 1) % ncount;
  double lhs = 0.0;
  int rhs = -1;

  for(int i = 1; i < ncount - 1; i++){
    int body_end = (root + i) % ncount;
    int new_vx = best_tour_nodes[body_end];
    SimpleTooth cand(root, body_start, body_end);

    increment_slack(cand, new_vx, lhs, rhs);

    //NOTE: we are completely ignoring heavy teeth for now
    if(cand.slack >= 0.5 - LP::EPSILON || cand.slack < 0)
      continue;

    if(body_size(cand) > (ncount - 2) / 2)
      complement(cand);

    if(body_size(cand) == 1){
      if(cand.root > cand.body_start &&
	 !light_teeth[cand.body_start].empty()){
	bool found_dup = false;
	for(auto orig = light_teeth[cand.body_start].rbegin();
	    orig != light_teeth[cand.body_start].rend(); orig++){
	  if(body_size(**orig) > 1) break;

	  if((*orig)->body_start == cand.root &&
	     cand.body_start == (*orig)->root){
	    found_dup = true;
	    break;
	  }
	}
	if(found_dup)
	  continue;
      }
    }

    try {
      light_teeth[root].emplace_back(SimpleTooth::Ptr(new SimpleTooth(cand)));
    } catch(...){
      PSEP_SET_GOTO(rval, "Couldn't push back light tooth. ");
    }
  }

  for(int k = 0; k < edge_marks.size(); k++) edge_marks[k] = 0;

 CLEANUP:
  if(rval)
    cerr << "CandidateTeeth::get_adjacent_teeth failed\n";
  return rval;
}

int CandidateTeeth::get_distant_teeth(const int root)
{
  int rval = 0;
  int ncount = G_s.node_count;

  for(int i = 2; i < ncount - 1; i++){
    int body_start = (root + i) % ncount;
    double lhs = 0.0;
    int rhs = -1;

    for(int j = i; j < ncount - 1; j++){
      int body_end = (root + j) % ncount;
      SimpleTooth cand(root, body_start, body_end);
      int new_vx = best_tour_nodes[body_end];

      increment_slack(cand, new_vx, lhs, rhs);

      if(cand.slack >= 0.5 - LP::EPSILON || cand.slack < 0)
	continue;

      if(body_size(cand) > (ncount - 2) / 2)
	complement(cand);

      if(body_size(cand) == 1){
	if(cand.root > cand.body_start &&
	   !light_teeth[cand.body_start].empty()){
	  bool found_dup = false;
	  for(auto orig = light_teeth[cand.body_start].rbegin();
	      orig != light_teeth[cand.body_start].rend(); orig++){
	    if(body_size(**orig) > 1) break;

	    if((*orig)->body_start == cand.root &&
	       cand.body_start == (*orig)->root){
	      found_dup = true;
	      break;
	    }
	  }
	  if(found_dup)
	    continue;
	}
      }

      try {
	light_teeth[root].emplace_back(SimpleTooth::Ptr(new SimpleTooth(cand)));
      } catch(...){
	PSEP_SET_GOTO(rval, "Couldn't push back light tooth. ");
      }
			    
    }
    for(int k = 0; k < edge_marks.size(); k++) edge_marks[k] = 0;
  }

 CLEANUP:
  if(rval)
    cerr << "CandidateTeeth::get_distant_teeth failed\n";
  return rval;
}

inline void CandidateTeeth::clear_collection()
{
  for(vector<SimpleTooth::Ptr> &vec : light_teeth){
#ifndef PSEP_TOOTH_UNIQ
    for(SimpleTooth::Ptr &T : vec)
      delete(T);
#endif
    vec.clear();
  }
}

inline bool SimpleTooth::sandwich() const
{
  return
    (body_start <= body_end) ?
    (body_start <= root && root <= body_end) : //[___<----*--->__]
    (body_start <= root || root <= body_end); // [-->__<--*--] OR [-*->__<--]
}

inline int CandidateTeeth::body_size(const SimpleTooth &T)
{
  return
    (T.body_start <= T.body_end) ?
    (T.body_end - T.body_start + !T.sandwich()) :
    ((G_s.node_count - T.body_start) + T.body_end + !T.sandwich());    
}

void CandidateTeeth::complement(SimpleTooth &T)
{
  int ncount = G_s.node_count;
  int body_start = T.body_start, body_end = T.body_end, root = T.root;
  int c_start, c_end;

  if((body_start == ((root + 1) % ncount)) ||
     (root == ((body_end + 1) % ncount))){
    if(body_start == ((root + 1) % ncount)){
      c_start = (body_end + 1) % ncount;
      c_end = (root + ncount - 1) % ncount;
    } else {
      c_start = (root + 1) % ncount;
      c_end = (body_end + ncount - 1) % ncount;
    }
  } else {
    c_start = (body_end + 1) % ncount;
    c_end = (body_start + ncount - 1) % ncount;
  }

  T.body_start = c_start;
  T.body_end = c_end;
}

int CandidateTeeth::body_subset(const SimpleTooth &T, const SimpleTooth &R,
				bool &result)
{
  if(T.root != R.root){
    cerr << "Cannot currently test body subset with different roots\n";
    return 1;
  }

  if(R.body_start <= R.body_end){ //if R is contiguous
    if(T.body_start > T.body_end){//R cannot contain a sandwich tooth T
      result = false;
    } else {
      result = (R.body_start <= T.body_start && T.body_end <= R.body_end);
    }
  } else { //if R wraps
    result = ( (R.body_start <= T.body_start &&
		(R.body_start <= T.body_end || T.body_end <= R.body_end)) ||
	       (T.body_start <= R.body_end && T.body_end <= R.body_end));
  }

  return 0;
}

//TODO/NOTES: it may be possible to speed this up somewhat with the edge
//indices hash map rather than the support graph
void CandidateTeeth::increment_slack(SimpleTooth &T, const int new_vx,
				     double &lhs, int &rhs)
{
  SNode *current_node = &(G_s.nodelist[new_vx]);

  rhs += 2;
  edge_marks[new_vx] = 1;

  for(int i = 0; i < current_node->s_degree; i++){
    int other_end = current_node->adj_objs[i].other_end;
    double lp_weight = current_node->adj_objs[i].lp_weight;

    if(edge_marks[other_end] == 1){
      lhs += (2 * lp_weight);
      continue;
    }

    if(other_end == best_tour_nodes[T.root]){
      lhs += lp_weight;
    }
  }

  T.slack = rhs - lhs;
}

void CandidateTeeth::print_tooth(const SimpleTooth &T)
{
  int ncount = G_s.node_count;
  int current_node = T.body_start;
  int upper_limit = body_size(T) + T.sandwich();

  cout << "Root: " << best_tour_nodes[T.root] << ", body size: "
       << body_size(T) 
       << ", Body: \n"
       << best_tour_nodes[T.body_start] << "\n";
  for(int i = 1; i < upper_limit; i++){
    current_node = best_tour_nodes[(T.body_start + i) % ncount];
    if(current_node != best_tour_nodes[T.root])
      cout <<current_node << "\n";
  }
  cout << "Slack: " << T.slack << "\n";
}

void CandidateTeeth::print_collection()
{
  for(int i = 0; i < G_s.node_count; i++){
    cout << "=== LIGHT TEETH WITH ROOT "
	 << best_tour_nodes[i] << ", " << light_teeth[i].size() << " TOTAL==="
	 << endl;
    // for(SimpleTooth::Ptr &T : light_teeth[i])
    //   print_tooth(*T);
  }
}
