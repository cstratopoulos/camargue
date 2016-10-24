#include "tooth.hpp"

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <iostream>
#include <algorithm>

using namespace PSEP;
using namespace std;

static int num_adjacent = 0, num_distant = 0;

CandidateTeeth::CandidateTeeth(vector<int> &_delta, vector<int> &_edge_marks,
			       vector<int> &_best_tour_nodes,
			       vector<int> &_perm,
			       SupportGraph &_G_s,
			       vector<int> &_support_indices,
			       vector<int> &_support_elist,
			       vector<double> &_support_ecap) :
  edge_marks(_edge_marks),
  best_tour_nodes(_best_tour_nodes),
  G_s(_G_s), support_elist(_support_elist), support_ecap(_support_ecap),
  cb_data(light_teeth,
	  _edge_marks,
	  _best_tour_nodes, _perm,
	  _G_s, _support_indices, _support_elist, _support_ecap)
{
  light_teeth.resize(best_tour_nodes.size());
}


int CandidateTeeth::get_light_teeth()
{
  int rval = 0;
  int notsort = 0;
  seg lin_seg(G_s.node_count - 1, G_s.node_count - 1, 0);
  vector<int> endmark;
  double ft, st;
  int max_deg = 0;

  
  try { endmark.resize(G_s.node_count, CC_LINSUB_BOTH_END); } catch(...) {
    PSEP_SET_GOTO(rval, "Couldn't set endmark. ");
  }

  for(int i = 0; i < G_s.node_count; i++)
    if(G_s.nodelist[i].s_degree > max_deg)
      max_deg = G_s.nodelist[i].s_degree;

  cout << "Support graph has max degree " << max_deg << "\n";

  cb_data.old_seg = &lin_seg;
  cb_data.cb_edge_marks[best_tour_nodes[lin_seg.start]] = 1;
  cb_data.unsorted_roots.clear();
  clear_collection();
  
  cout << "Getting light teeth via linsub...." << endl;
  ft = zeit();
  
  rval = CCcut_linsub_allcuts(G_s.node_count,
			      G_s.edge_count,
			      &best_tour_nodes[0], &endmark[0],
			      &support_elist[0], &support_ecap[0], 2.999,
			      &cb_data, get_teeth);
  if(rval) goto CLEANUP;
  ft = zeit() - ft;

  cout << ft << "s total, " << num_adjacent << " adjacent teeth, "
       << num_distant << " distant teeth\n";

  st = zeit();
  for(int i = 0; i < cb_data.unsorted_roots.size(); i++){
    int root = cb_data.unsorted_roots[i];
    std::sort(light_teeth[root].begin(), light_teeth[root].end(),
	      [this](const SimpleTooth::Ptr &T,
		     const SimpleTooth::Ptr &R) -> bool {
		return body_size(*T) < body_size(*R);
	      });
  }
  st = zeit() - st;

  cout << "Sorted only " << cb_data.unsorted_roots.size() << " lists in "
       << st << "s\n";


  cout << "Checking sorted status....";
  for(vector<SimpleTooth::Ptr> &t_vec : light_teeth){
    bool is_s = std::is_sorted(t_vec.begin(), t_vec.end(),
			       [this](const SimpleTooth::Ptr &T,
				      const SimpleTooth::Ptr &R) -> bool {
				 return body_size(*T) < body_size(*R);
			       });
    if(!is_s){
      notsort++;
      cout << "List of teeth with root "
	   << t_vec.front()->root << " is not sorted\n";
    }
  }

  if(notsort == 0)
    cout << "All vectors sorted by increasing body size\n";
  else
    cout << notsort << " vectors not sorted\n";

 CLEANUP:
  if(rval == 1){
    cerr << "CandidateTeeth::get_light_teeth failed\n";
    clear_collection();
  }
  return rval;
}

int CandidateTeeth::get_teeth(double cut_val, int cut_start, int cut_end,
			      void *u_data)
{
  int rval = 0;

  LinsubCBData *arg = (LinsubCBData *) u_data;
  
  seg *old_cut = arg->old_seg;
  vector<vector<SimpleTooth::Ptr>> &teeth = arg->cb_teeth;
  
  vector<int> &marks = arg->cb_edge_marks;

  vector<int> &best_nodes = arg->cb_tour_nodes;
  vector<int> &perm = arg->cb_perm;

  SupportGraph &G = arg->cb_G_s;

  unordered_map<int, double> &rb_sums = arg->root_bod_sums;

  int ncount = best_nodes.size();
  double slack = (cut_val - 2.0) / 2.0;
  int set_size = cut_end - cut_start + 1;
  int rhs = (2 * set_size) - 1;
  double partial_lhs = (2 * set_size) - cut_val;
  double root_bod_lb = rhs - partial_lhs - 0.4999;

  if(cut_start == old_cut->start){//if the current seg contains previous
    if(cut_end == old_cut->end + 1 && slack + old_cut->cutval < 0.4999){
      rval = add_tooth(teeth, cut_end, cut_start, old_cut->end,
		       old_cut->cutval + slack);
      PSEP_CHECK_RVAL(rval, "Problem with adjacent teeth. ");
      num_adjacent++;
    }
    
    for(int i = old_cut->end + 1; i <= cut_end; i++){
      marks[best_nodes[i]] = 1;
      rb_sums.erase(i);
    }
    rb_sums.erase((cut_end + 1) % ncount);
    marks[best_nodes[(cut_end + 1) % ncount]] = 1;
    
	  

    for(int i = old_cut->end + 1; i <= cut_end; i++){
      SNode current_node = G.nodelist[best_nodes[i]];

      for(int k = 0; k < current_node.s_degree; k++){
    	int other_end = current_node.adj_objs[k].other_end;
    	int root_perm = perm[other_end];

    	if(marks[other_end] == 0)//if not in or adjacent to body
    	  rb_sums[root_perm] += current_node.adj_objs[k].lp_weight;
      }
    }
    
  } else { //if we are on a new start (degree eqn "segment")
    for(int i = old_cut->end + 1; i <= cut_end; i++)
      marks[best_nodes[i]] = 0;
    marks[best_nodes[cut_start]] = 1;
    marks[best_nodes[(cut_start + 1) % ncount]] = 1;
    marks[best_nodes[(cut_start + ncount - 1) % ncount]] = 1;

    rb_sums.clear();

    SNode current_node = G.nodelist[best_nodes[cut_start]];
    
    for(int k = 0; k < current_node.s_degree; k++){
      int other_end = current_node.adj_objs[k].other_end;
      int root_perm = perm[other_end];

      if(marks[other_end] == 0)//if not in or adjacent to body
  	rb_sums[root_perm] += current_node.adj_objs[k].lp_weight;
    }
  }

  for(auto &kv : rb_sums)
    if(kv.second > root_bod_lb){
      rval = add_tooth(teeth, kv.first, cut_start, cut_end,
  		       rhs - partial_lhs - kv.second);
      PSEP_CHECK_RVAL(rval, "Problem with distant teeth. ");
      
      if(arg->unsorted_roots.empty() || arg->unsorted_roots.back() != kv.first)
	arg->unsorted_roots.push_back(kv.first);
      num_distant++;
    }
  


 CLEANUP:
  if(rval)
    cerr << "Linsub callback failed\n";
  old_cut->start = cut_start;
  old_cut->end = cut_end;
  old_cut->cutval = slack;
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

inline int CandidateTeeth::add_tooth(vector<vector<SimpleTooth::Ptr>> &teeth,
				     const int root, const int body_start,
				     const int body_end, const double slack)
{
  int rval = 0;

  try {
    teeth[root].emplace_back(SimpleTooth::Ptr(new SimpleTooth(root,
							      body_start,
							      body_end,
							      slack)));
  } catch(...) { rval = 1; }

  if(rval)
    cerr << "Couldn't push back light tooth\n";
  return rval;
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

inline void CandidateTeeth::LinsubCBData::refresh(seg *new_old_seg)
{
  new_old_seg->start = cb_tour_nodes.size() - 1;
  new_old_seg->end = new_old_seg->start;
  new_old_seg->cutval = 0;
  old_seg = new_old_seg;
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
