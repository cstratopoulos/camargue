#include "n_tooth.hpp"

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <iostream>
#include <algorithm>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::pair;
using std::unique_ptr;
using std::inplace_merge;

namespace PSEP {
namespace nu {

static bool ptr_cmp(const SimpleTooth::Ptr &S, const SimpleTooth::Ptr &T)
{ return S->body_size() < T->body_size(); }

CandidateTeeth::CandidateTeeth(PSEP::Data::GraphGroup &_graph_dat,
			       PSEP::Data::BestGroup &_best_dat,
			       PSEP::Data::SupportGroup &_supp_dat) :
  light_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  left_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  right_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  dist_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  adj_zones(vector<vector<int>>(_supp_dat.G_s.node_count,
				vector<int>(_supp_dat.G_s.node_count, 0))),
  stats(_supp_dat.G_s.node_count, ListStat::None),
  endmark(vector<int>(_supp_dat.G_s.node_count, CC_LINSUB_BOTH_END)),
  graph_dat(_graph_dat),
  best_dat(_best_dat),
  supp_dat(_supp_dat)
{
  PSEP::SupportGraph &G_s = supp_dat.G_s;
  int ncount = G_s.node_count;
  vector<int> &perm = best_dat.perm;
  vector<int> &tour = best_dat.best_tour_nodes;

  #pragma omp parallel for
  for(int root_ind = 0; root_ind < ncount; ++root_ind){
    int actual_vx = tour[root_ind];
    PSEP::SNode x = G_s.nodelist[actual_vx];
    
    for(int k = 0; k < x.s_degree; ++k){
      int end_ind = perm[x.adj_objs[k].other_end];
      
      adj_zones[root_ind][end_ind] = 1;
    }

    int label = 0;
    for(int &i : adj_zones[root_ind]){
      if(i == 1){
	++label;
	i = -1 * label;
      } else
	i = label;
    }
  }
}

int CandidateTeeth::get_light_teeth()
{
  int rval = 0;
  unique_ptr<LinsubCBData> cb_data;

  try {
    cb_data =
      PSEP::make_unique<LinsubCBData>(right_teeth,left_teeth, dist_teeth,
				      adj_zones, graph_dat.edge_marks,
				      best_dat.best_tour_nodes,
				      best_dat.perm, supp_dat.G_s); }
  catch(...){ PSEP_SET_GOTO(rval, "Couldn't allocate CBdata. "); }

  rval = CCcut_linsub_allcuts(supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
			      &best_dat.best_tour_nodes[0], &endmark[0],
			      &supp_dat.support_elist[0],
			      &supp_dat.support_ecap[0], 3.0 - Epsilon::Cut,
			      cb_data.get(), teeth_cb);
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval)
    cerr << "Problem in CandidateTeeth::get_light_teeth.\n";
  return rval;
}

int CandidateTeeth::merge_and_sort(const int root)
{
  vector<SimpleTooth::Ptr> &teeth = light_teeth[root];
  int left_sz = left_teeth[root].size();
  int right_sz = right_teeth[root].size();
  int dist_sz = dist_teeth[root].size();
  
  try {
    for(SimpleTooth::Ptr &T : left_teeth[root]) teeth.push_back(std::move(T));
    for(SimpleTooth::Ptr &T : right_teeth[root]) teeth.push_back(std::move(T));
    for(SimpleTooth::Ptr &T : dist_teeth[root]) teeth.push_back(std::move(T));
  } catch(...){ cerr << "CandidateTeeth::merge_and sort failed.\n"; return 1; }

  if(dist_sz > 0){
    std::sort(teeth.begin(), teeth.end(), ptr_cmp);
    stats[root] = ListStat::Full;
  } else
    if(left_sz > 0 && right_sz > 0){
      std::inplace_merge(teeth.begin(), teeth.begin() + left_sz, teeth.end(),
			 ptr_cmp);
      stats[root] = ListStat::Merge;
    }
  return 0;
}

void CandidateTeeth::weak_elim()
{
  for(int root = 0; root < light_teeth.size(); ++root){
    if(stats[root] == ListStat::None) continue;
    
    vector<SimpleTooth::Ptr> &teeth = light_teeth[root];
    bool found_elim = false;
    
    for(auto it = teeth.begin(); it != teeth.end() - 1; ++it){
      SimpleTooth::Ptr &S = *it;      
      if(S->root == -1) continue;
      
      for(auto it2 = it + 1; it2 != teeth.end(); ++it2){
	SimpleTooth::Ptr &T = *it2;	
	if(T->root == -1) continue;

	if(root_equivalent(root, tooth_seg(S->body_start, S->body_end),
			   tooth_seg(T->body_start, T->body_end))){
	  found_elim = true;
	  if(S->slack < T->slack)
	    T->root = -1;
	  else
	    S->root = -1;
	}
      }
    }
    
    if(found_elim)
      teeth.erase(std::remove_if(teeth.begin(), teeth.end(),
				 [](const SimpleTooth::Ptr &T) -> bool {
				   return T->root == -1;
				 }),
		  teeth.end());
  }
}

void CandidateTeeth::get_range(const int root, const tooth_seg &s,
			       IntPair &range,
			       const vector<vector<int>> &zones)
{
  range = IntPair(-1, -1);
  int start = zones[root][s.start], end = zones[root][s.end];

  if(start == end){//if same endpoints
    if(start < 0){//singleton
      range = IntPair(-start, -start);
    }//else empty
    return;
  }

  //now different endpoints
  if(start < 0){
    range.first = -start; //first definitely start
    range.second = fabs(end);
    return;
  }

  //now diff endpoints, start >= 0
  if(end < 0){
    range.second = -end; //second definitely end
    if(start == -end)
      range.first = start; //start is itself
    else
      range.first = start + 1;// fmin(start + 1, -end) start is at most end
    return;
  }

  //now diff endpoints, both >= 0
  range = IntPair(start + 1, end);// fmax(start + 1, end - 1));
}

bool CandidateTeeth::root_equivalent(const int root, const tooth_seg &s1,
				     const tooth_seg &s2,
				     const vector<vector<int>> &zones)
{
  IntPair s1_range, s2_range;
    
  get_range(root, s1, s1_range, zones);
  get_range(root, s2, s2_range, zones);
  return s1_range == s2_range;
}

bool CandidateTeeth::root_equivalent(const int root, const tooth_seg &s1,
				     const tooth_seg &s2) const
{
  return root_equivalent(root, s1, s2, adj_zones);
}

int CandidateTeeth::teeth_cb(double cut_val, int cut_start, int cut_end,
			     void *u_data)
{
  int rval = 0;
  double slack = (cut_val - 2.0) / 2.0;
  
  LinsubCBData *arg = (LinsubCBData *) u_data;
  
  //distant declarations
  vector<int>
    &marks = arg->node_marks, &tour = arg->tour_nodes, &perm = arg->perm;
  PSEP::SupportGraph &G = arg->G_s;
  std::unordered_map<int, double> &rb_sums = arg->rb_sums;
  int ncount = G.node_count;
  double rb_lower = cut_val - (1.5 - Epsilon::Cut);

  //right-adjacent declarations
  tooth_seg &old_seg = arg->old_seg;
  vector<vector<int>> &zones = arg->adj_zones;
  vector<pair<int, double>> &old_rights = arg->prev_slacks;

  //right adjacent add/elim
  if(cut_start == old_seg.start){
    if(cut_end == old_seg.end + 1 &&
       slack + old_seg.slack < (0.5 - Epsilon::Cut)){
      int root = cut_end;
      double new_slack = slack + old_seg.slack;
      vector<SimpleTooth::Ptr> &r_vec = arg->r_teeth[root];
      bool elim = false;
      
      if(!r_vec.empty()){
	tooth_seg prev_body(r_vec.back()->body_start, r_vec.back()->body_end);
	double prev_slack = r_vec.back()->slack;
	if(CandidateTeeth::root_equivalent(root, prev_body, old_seg, zones)){
	  elim = true;
	  if(new_slack < prev_slack)
	    r_vec.back() = PSEP::make_unique<SimpleTooth>(root, old_seg,
							  new_slack);
	}
      }

      if(!elim){
	try {
	  r_vec.emplace_back(PSEP::make_unique<SimpleTooth>(root, old_seg,
							    new_slack));
	} catch(...){ PSEP_SET_GOTO(rval, "Couldn't push back new tooth. ")}
      }
    }
  }

  //left adjacent add/elim
  if(cut_start + 1 != cut_end){
    pair<int, double> &old_right_pair = arg->prev_slacks[cut_end];
    
    if((old_right_pair.first == cut_start + 1) &&
       (old_right_pair.second + slack < (0.5 - Epsilon::Cut))){
      int root = cut_start;
      double new_slack = slack + old_right_pair.second;
      tooth_seg new_body(old_right_pair.first, cut_end);
      vector<SimpleTooth::Ptr> &l_vec = arg->l_teeth[root];
      bool elim = false;

      if(!l_vec.empty()){
	tooth_seg prev_body(l_vec.back()->body_start, l_vec.back()->body_end);
	double prev_slack = l_vec.back()->slack;

	if(CandidateTeeth::root_equivalent(root, new_body, prev_body, zones)){
	  elim = true;
	  if(new_slack < prev_slack)
	    l_vec.back() = PSEP::make_unique<SimpleTooth>(root, new_body,
							  new_slack);
	}
      }

      if(!elim){
	try {
	  l_vec.emplace_back(PSEP::make_unique<SimpleTooth>(root, new_body,
							    new_slack));
	} catch(...){ PSEP_SET_GOTO(rval, "Couldn't push back new tooth. "); }
      }
    }
  }
  

  //distant add
  if(cut_start == old_seg.start){
    for(int i = old_seg.end + 1; i <= cut_end; ++i){
      marks[i] = 1;
      rb_sums.erase(i);
    }
    marks[(cut_end + 1) % ncount] = 1;
    rb_sums.erase((cut_end + 1) % ncount);

    for(int i = old_seg.end + 1; i <= cut_end; ++i){
      SNode vx = G.nodelist[tour[i]];
      for(int k = 0; k < vx.s_degree; ++k){
	int root_perm = perm[vx.adj_objs[k].other_end];

	if(marks[root_perm] == 0)
	  rb_sums[root_perm] += vx.adj_objs[k].lp_weight;
      }
    }
  } else {//on a new degree node
    //clean up after the old one
    for(int i = 0; i < ncount; ++i){
      marks[i] = 0;
    }
    rb_sums.clear();

    //set up for the new one
    marks[cut_start] = 1;
    marks[(cut_start + 1) % ncount] = 1;
    marks[(cut_start + (ncount - 1)) % ncount] = 1;

    SNode vx = G.nodelist[tour[cut_start]];
    for(int k = 0; k < vx.s_degree; ++k){
      int root_perm = perm[vx.adj_objs[k].other_end];

      if(marks[root_perm] == 0)
	rb_sums[root_perm] += vx.adj_objs[k].lp_weight;
    }    
  }

  for(auto &kv : rb_sums){
    int i = kv.first;
    double rb_sum = kv.second;
    if(rb_sum > rb_lower){
      vector<SimpleTooth::Ptr> &dt = arg->d_teeth[i];
      bool elim = false;
      double new_slack = cut_val - rb_sum - 1;
      tooth_seg new_body(cut_start, cut_end);

      if(!dt.empty()){
	tooth_seg prev_body(dt.back()->body_start, dt.back()->body_end);
	double prev_slack = dt.back()->slack;

	if(CandidateTeeth::root_equivalent(i, new_body, prev_body, zones)){
	  elim = true;
	  if(new_slack < prev_slack)
	    dt.back() = PSEP::make_unique<SimpleTooth>(i, new_body, new_slack);
	}
      }
      
      if(!elim){
	try {
	  dt.emplace_back(PSEP::make_unique<SimpleTooth>(i, new_body,
							 new_slack));
	} catch(...){ PSEP_SET_GOTO(rval, "Couldn't push back tooth. "); }
      }
    }
  }


 CLEANUP:
  if(rval)
    cerr << "CandidateTeeth::tooth_cb failed.\n";
  old_seg = tooth_seg(cut_start, cut_end, slack); //right adjacent update
  old_rights[cut_end] = pair<int, double>(cut_start, slack); //left adj update
  return rval;
}

string CandidateTeeth::print_label(const SimpleTooth &T)
{
  string
    body = (T.body_start == T.body_end) ?
    to_string(T.body_start) :
    ("{"
     + to_string(T.body_start) + "..." + to_string(T.body_end)
     + "}");

  return "(" + to_string(T.root) + ", " + body + ")";
}

void CandidateTeeth::print_tooth(const SimpleTooth &T, bool full)
{
  print_tooth(T, full, best_dat.best_tour_nodes);
}

void CandidateTeeth::print_tooth(const SimpleTooth &T, bool full,
				 const vector<int> &bt)
{
  cout << "(" << bt[T.root] << ", {" << bt[T.body_start];

  if(T.body_start == T.body_end){
    cout << "}) -- slack " << T.slack << "\n";
    return;
  }

  if(!full){
    cout << ", ..., " << bt[T.body_end] << "}) -- slack " << T.slack << "\n";
    return;
  }

  bool comma_sep = T.body_end - T.body_start <= 20;
  int i = T.body_start;
  
  cout << ", ";
  while(i++ != T.body_end){
    cout << bt[i];
    if(i == T.body_end)
      cout << "}";
    else
      cout << (comma_sep ? ", " : "\n\t");
  }
  cout << ") -- slack " << T.slack << "\n";  
}

}
}
