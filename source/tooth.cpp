#include "tooth.hpp"

#include <timsort.hpp>

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <iostream>
#include <algorithm>
#include <tuple>

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

static inline bool ptr_cmp(const SimpleTooth::Ptr &S, const SimpleTooth::Ptr &T)
{ return S->body_size() < T->body_size(); }

static inline bool ptr_elim(const SimpleTooth::Ptr &S){ return S->root == -1; }

static inline bool elim_less_tie(const SimpleTooth::Ptr &S,
				 const SimpleTooth::Ptr &T)
{
  return std::make_tuple(S->slack, S->body_size()) <
    std::make_tuple(T->slack, T->body_size());
}

CandidateTeeth::CandidateTeeth(Data::GraphGroup &_graph_dat,
			       Data::BestGroup &_best_dat,
			       Data::SupportGroup &_supp_dat) :
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
  supp_dat(_supp_dat),
  t_all("Candidate Teeth"),
  t_zones("Adj zones", &t_all),
  t_find("Initial find", &t_all),
  t_elim("Weak elim", &t_all),
  t_sort("Merge and sort", &t_all)
{
  t_all.start();
  t_zones.start();
  SupportGraph &G_s = supp_dat.G_s;
  int ncount = G_s.node_count;
  vector<int> &perm = best_dat.perm;
  vector<int> &tour = best_dat.best_tour_nodes;

  #pragma omp parallel for
  for(int root_ind = 0; root_ind < ncount; ++root_ind){
    int actual_vx = tour[root_ind];
    SNode x = G_s.nodelist[actual_vx];
    
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
  
  t_zones.stop();
}

int CandidateTeeth::get_light_teeth()
{
  t_find.start();
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
  t_find.stop();

 CLEANUP:
  if(rval)
    cerr << "Problem in CandidateTeeth::get_light_teeth.\n";
  return rval;
}

int CandidateTeeth::merge_and_sort()
{
  t_sort.start();
  int rval = 0;

  #pragma omp parallel for
  for(int root = 0; root < light_teeth.size(); ++root){
    if(rval) continue;
    if(merge_and_sort(root)){
      #pragma omp critical
      { rval = 1; }
    }
  }
  
  t_sort.stop();
  t_all.stop();
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
  } catch(...){
    cerr << "CandidateTeeth::merge_and sort failed.\n"; return 1;
  }

  if(dist_sz > 0 || (left_sz > 0 && right_sz > 0)){
    gfx::timsort(teeth.begin(), teeth.end(), ptr_cmp);
    stats[root] = ListStat::Full;
  }

  return 0;
}

void CandidateTeeth::unmerged_weak_elim()
{
  t_elim.start();

  #pragma omp parallel for
  for(int root = 0; root < light_teeth.size(); ++root){
    vector<SimpleTooth::Ptr>
      &right = right_teeth[root],
      &left = left_teeth[root],
      &dist = dist_teeth[root];

    bool right_elim = false;
    bool left_elim = false;
    bool dist_elim = false;

    if(dist.empty()) continue;

    for(SimpleTooth::Ptr &D : dist){
      if(D->root == -1) continue;
      
      if(D->root < D->body_start){//left of body
	for(SimpleTooth::Ptr &L : left){
	  if(L->root == -1) continue;
	  if(D->root == -1) break;
	  if(root_equivalent(root, tooth_seg(L->body_start, L->body_end),
			     tooth_seg(D->body_start, D->body_end))){
	    if(elim_less_tie(D, L)){
	      L->root = -1;
	      left_elim = true;
	    } else {
	      D->root = -1;
	      dist_elim = true;
	    }
	  }
	}
	continue;
      }

      if(D->root > D->body_end){//right of body
	for(SimpleTooth::Ptr &R : right){
	  if(R->root == -1) continue;
	  if(D->root == -1) break;
	  if(root_equivalent(root, tooth_seg(R->body_start, R->body_end),
			     tooth_seg(D->body_start, D->body_end))){
	    if(elim_less_tie(D,R)){
	      R->root = -1;
	      right_elim = true;
	    } else {
	      D->root = -1;
	      dist_elim = true;
	    }
	  }	  
	}
      }
    }

    if(right_elim)
      right.erase(std::remove_if(right.begin(), right.end(), ptr_elim),
		  right.end());
    if(left_elim)
      left.erase(std::remove_if(left.begin(), left.end(), ptr_elim),
		 left.end());
    if(dist_elim)
      dist.erase(std::remove_if(dist.begin(), dist.end(), ptr_elim),
		 dist.end());    
  }
  
  t_elim.stop();
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

  //discard negative values caused by numerical instability
  double slack = (cut_val - (2.0 - Epsilon::Cut)) / 2.0;
  
  LinsubCBData *arg = (LinsubCBData *) u_data;
  if((cut_end - cut_start + 1) > (arg->G_s.node_count - 2)) return 0;
  
  //distant declarations
  vector<int>
    &marks = arg->node_marks, &tour = arg->tour_nodes, &perm = arg->perm;
  SupportGraph &G = arg->G_s;
  std::unordered_map<int, double> &rb_sums = arg->rb_sums;
  int ncount = G.node_count;
  double rb_lower = cut_val - (1.5 - Epsilon::Cut);

  //right-adjacent declarations
  tooth_seg &old_seg = arg->old_seg;
  vector<vector<int>> &zones = arg->adj_zones;
  vector<pair<int, double>> &old_rights = arg->prev_slacks;

  //distant add
  if(cut_start == old_seg.start){
    for(int i = old_seg.end + 1; i <= cut_end; ++i){
      marks[i] = 1;
      rb_sums.erase(i);
    }

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

    SNode vx = G.nodelist[tour[cut_start]];
    for(int k = 0; k < vx.s_degree; ++k){
      int root_perm = perm[vx.adj_objs[k].other_end];

      if(marks[root_perm] == 0)
	rb_sums[root_perm] += vx.adj_objs[k].lp_weight;
    }    
  }

  for(auto &kv : rb_sums){
    int root = kv.first;
    double rb_sum = kv.second;
    if((cut_start == cut_end) && (root > cut_end) && (root != (ncount - 1)))
      continue;
    if(((cut_end - cut_start + 1) == (ncount - 2)) && (root > cut_end))
      continue;

    if(rb_sum > rb_lower){
      vector<SimpleTooth::Ptr>
	&teeth = (root + 1 == cut_start) ? arg->l_teeth[root] :
	((cut_end + 1 == root) ? arg->r_teeth[root] : arg->d_teeth[root]);
      double abs_slack = fabs(cut_val - rb_sum - 1);
      double new_slack = (abs_slack < Epsilon::Zero) ? 0 : abs_slack;

      try {
	add_tooth(teeth, zones, root, cut_start, cut_end, new_slack);
      } catch(...){ PSEP_SET_GOTO(rval, "Couldn't push back dist tooth. "); }
    }
  }


 CLEANUP:
  if(rval)
    cerr << "CandidateTeeth::tooth_cb failed.\n";
  old_seg = tooth_seg(cut_start, cut_end, slack); //right adjacent update
  old_rights[cut_end] = pair<int, double>(cut_start, slack); //left adj update
  return rval;
}

inline void CandidateTeeth::add_tooth(vector<SimpleTooth::Ptr> &teeth,
			      const vector<vector<int>> &zones,
			      const int root, const int body_start,
			      const int body_end, const double slack)
{
  bool elim = false;
  tooth_seg body(body_start, body_end);
  
  if(!teeth.empty()){
    tooth_seg old_body(teeth.back()->body_start, teeth.back()->body_end);
    double old_slack{teeth.back()->slack};
    if(CandidateTeeth::root_equivalent(root, body, old_body, zones)){
      elim = true;
      if(slack < old_slack)
	teeth.back() = PSEP::make_unique<SimpleTooth>(root, body, slack);
    }
  }

  if(!elim)
    teeth.emplace_back(PSEP::make_unique<SimpleTooth>(root, body, slack));
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

void CandidateTeeth::profile()
{
  t_zones.report(true);
  t_find.report(false);
  t_elim.report(true);
  t_sort.report(true);
  t_all.report(true);
}

}
