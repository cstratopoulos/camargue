#include "tooth.hpp"
#include "config.hpp"

#if CMR_HAVE_TIMSORT
#include <timsort.hpp>
#endif

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

namespace CMR {

using ToothList = CandidateTeeth::ToothList;

static inline bool ptr_cmp(const SimpleTooth::Ptr &S, const SimpleTooth::Ptr &T)
{ return S->body_size() < T->body_size(); }

static inline bool ptr_elim(const SimpleTooth::Ptr &S) { return S->root == -1; }

static inline bool elim_less_tie(const SimpleTooth::Ptr &S,
				 const SimpleTooth::Ptr &T)
{
  return std::make_tuple(S->slack, S->body_size()) <
    std::make_tuple(T->slack, T->body_size());
}

static void tooth_sort(ToothList &T)
{
#if CMR_HAVE_TIMSORT
    gfx::timsort(T.begin(), T.end(), ptr_cmp);
#else
    std::sort(T.begin(), T.end(), ptr_cmp);
#endif
}

vector<vector<int>> CandidateTeeth::adj_zones;
vector<util::SquareUT<ToothList::reverse_iterator>> CandidateTeeth::seen_ranges;

CandidateTeeth::CandidateTeeth(Data::GraphGroup &_graph_dat,
			       Data::BestGroup &_best_dat,
			       Data::SupportGroup &_supp_dat) :
  light_teeth(std::vector<ToothList>(_supp_dat.G_s.node_count)),
  left_teeth(std::vector<ToothList>(_supp_dat.G_s.node_count)),
  right_teeth(std::vector<ToothList>(_supp_dat.G_s.node_count)),
  dist_teeth(std::vector<ToothList>(_supp_dat.G_s.node_count)),
  stats(_supp_dat.G_s.node_count, ListStat::None),
  endmark(vector<int>(_supp_dat.G_s.node_count, CC_LINSUB_BOTH_END)),
  graph_dat(_graph_dat),
  best_dat(_best_dat),
  supp_dat(_supp_dat),
  t_all("Candidate Teeth"),
  t_zones("Adj zones", &t_all),
  t_find("Initial find", &t_all),
  t_elim("Weak elim", &t_all),
  t_comp_elim("Complement elim", &t_all),
  t_sort("Merge and sort", &t_all)
{
  t_all.start();
  t_zones.start();
  SupportGraph &G_s = supp_dat.G_s;
  int ncount = G_s.node_count;
  vector<int> &perm = best_dat.perm;
  vector<int> &tour = best_dat.best_tour_nodes;

  seen_ranges.resize(ncount);
  adj_zones.resize(ncount);
  for (vector<int> &vec : adj_zones) {
      vec.resize(ncount);
      for (int &i : vec)
          i = 0;
  }

#ifdef CMR_USE_OMP
  #pragma omp parallel for
#endif
  for (int root_ind = 0; root_ind < ncount; ++root_ind) {
    int actual_vx = tour[root_ind];
    SNode x = G_s.nodelist[actual_vx];
    light_teeth[root_ind].reserve(2 * (x.s_degree - 1));
    
    for (int k = 0; k < x.s_degree; ++k) {
      int end_ind = perm[x.adj_objs[k].other_end];
      
      adj_zones[root_ind][end_ind] = 1;
    }

    int label = 0;
    for (int &i : adj_zones[root_ind]) {
      if (i == 1) {
	++label;
	i = -1 * label;
      } else
	i = label;
    }

    seen_ranges[root_ind] =
    util::SquareUT<ToothList::reverse_iterator>(label + 1,
                                        light_teeth[root_ind].rend());
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
    util::make_unique<LinsubCBData>(light_teeth,
                                    adj_zones, seen_ranges,
                                    graph_dat.node_marks,
                                    best_dat.best_tour_nodes,
                                    best_dat.perm, supp_dat.G_s); }
  catch (...) { CMR_SET_GOTO(rval, "Couldn't allocate CBdata. "); }

  rval = CCcut_linsub_allcuts(supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
			      &best_dat.best_tour_nodes[0], &endmark[0],
			      &supp_dat.support_elist[0],
			      &supp_dat.support_ecap[0],
			      3.0 - Epsilon::Cut,
			      cb_data.get(), teeth_cb);
  if (rval) goto CLEANUP;
  t_find.stop();

 CLEANUP:
  if (rval)
    cerr << "Problem in CandidateTeeth::get_light_teeth.\n";
  return rval;
}

int CandidateTeeth::merge_and_sort()
{
  t_sort.start();
  int rval = 0;

#ifdef CMR_USE_OMP
  #pragma omp parallel for
#endif
  for (int root = 0; root < light_teeth.size(); ++root) {
    if (rval) continue;
    if (merge_and_sort(root)) {
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
    ToothList &teeth = light_teeth[root];
    tooth_sort(teeth);
    stats[root] = ListStat::Full;

  return 0;
}

void CandidateTeeth::complement_elim()
{
  t_comp_elim.start();
//   int ncount = light_teeth.size();
  
//   int root = ncount - 1;
//   if (!right_teeth[root].empty() && !dist_teeth[root].empty()) {
//     ToothList
//       &right = right_teeth[root], &dist = dist_teeth[root];
//     double right_elim = false;
//     double dist_elim = false;

//     for (SimpleTooth::Ptr &R : right) {
//       for (SimpleTooth::Ptr &D : dist) {
// 	if ((D->body_start == 0) && ((D->body_end + 1) == R->body_start)) {
// 	  if (R->body_size() <= D->body_size()) {
// 	    dist_elim = true;
// 	    D->root = -1;
// 	  } else {
// 	    right_elim = true;
// 	    R->root = -1;
// 	  }
// 	}
// 	if (right_elim) break;
//       }
//     }

//     if (right_elim)
//       right.erase(std::remove_if (right.begin(), right.end(), ptr_elim),
// 		  right.end());
//     if (dist_elim)
//       dist.erase(std::remove_if (dist.begin(), dist.end(), ptr_elim),
// 		 dist.end());
//   }

//   root = 0;
//   if (!left_teeth[root].empty() && !dist_teeth[root].empty()) {
//     ToothList
//       &left = left_teeth[root], &dist = dist_teeth[root];
//     bool left_elim = false;
//     bool dist_elim = false;

//     for (SimpleTooth::Ptr &L : left) {
//       for (SimpleTooth::Ptr &D : dist) {
// 	if ((D->body_end == (ncount - 1)) &&
// 	   ((L->body_end + 1)  == D->body_start)) {
// 	  if (L->body_size() <= D->body_size()) {
// 	    dist_elim = true;
// 	    D->root = -1;
// 	  } else {
// 	    left_elim = true;
// 	    L->root = -1;
// 	  }
// 	}
// 	if (left_elim) break;
//       }
//     }

//     if (left_elim)
//       left.erase(std::remove_if (left.begin(), left.end(), ptr_elim),
// 		 left.end());
//     if (dist_elim)
//       dist.erase(std::remove_if (dist.begin(), dist.end(), ptr_elim),
// 		 dist.end());
//   }

// #ifdef CMR_USE_OMP
//   #pragma omp parallel for
// #endif
//   for (root = 1; root < ncount - 1; ++root) {
//     ToothList
//       &right = right_teeth[root], &left = left_teeth[root];
//     if (right.empty() || left.empty()) continue;
//     bool right_elim = false;
//     bool left_elim = false;

//     for (SimpleTooth::Ptr &R : right) {
//       for (SimpleTooth::Ptr &L : left) {
// 	if ((L->body_end == (ncount - 1)) && (R->body_start == 0)) {
// 	  if (L->body_size() <= R->body_size()) {
// 	    right_elim = true;
// 	    R->root = -1;
// 	  } else {
// 	    left_elim = true;
// 	    L->root = -1;
// 	  }
// 	}
// 	if (right_elim) break;
//       }
//     }
    
//     if (left_elim)
//       left.erase(std::remove_if (left.begin(), left.end(), ptr_elim),
// 		 left.end());
//     if (right_elim)
//       right.erase(std::remove_if (right.begin(), right.end(), ptr_elim),
// 		  right.end());    
//   }
  
  t_comp_elim.stop();
}

void CandidateTeeth::unmerged_weak_elim()
{
  t_elim.start();
// #ifdef CMR_USE_OMP
//   #pragma omp parallel for
// #endif
//   for (int root = 0; root < light_teeth.size(); ++root) {
//     ToothList
//       &right = right_teeth[root],
//       &left = left_teeth[root],
//       &dist = dist_teeth[root];

//     bool right_elim = false;
//     bool left_elim = false;
//     bool dist_elim = false;

//     if (dist.empty()) continue;

//     for (SimpleTooth::Ptr &D : dist) {
//       if (D->root == -1) continue;
      
//       if (D->root < D->body_start) {//left of body
// 	for (SimpleTooth::Ptr &L : left) {
// 	  if (L->root == -1) continue;
// 	  if (D->root == -1) break;
// 	  if (root_equivalent(root, tooth_seg(L->body_start, L->body_end),
// 			     tooth_seg(D->body_start, D->body_end))) {
// 	    if (elim_less_tie(D, L)) {
// 	      L->root = -1;
// 	      left_elim = true;
// 	    } else {
// 	      D->root = -1;
// 	      dist_elim = true;
// 	    }
// 	  }
// 	}
// 	continue;
//       }

//       if (D->root > D->body_end) {//right of body
// 	for (SimpleTooth::Ptr &R : right) {
// 	  if (R->root == -1) continue;
// 	  if (D->root == -1) break;
// 	  if (root_equivalent(root, tooth_seg(R->body_start, R->body_end),
// 			     tooth_seg(D->body_start, D->body_end))) {
// 	    if (elim_less_tie(D,R)) {
// 	      R->root = -1;
// 	      right_elim = true;
// 	    } else {
// 	      D->root = -1;
// 	      dist_elim = true;
// 	    }
// 	  }	  
// 	}
//       }
//     }

//     if (right_elim)
//       right.erase(std::remove_if (right.begin(), right.end(), ptr_elim),
// 		  right.end());
//     if (left_elim)
//       left.erase(std::remove_if (left.begin(), left.end(), ptr_elim),
// 		 left.end());
//     if (dist_elim)
//       dist.erase(std::remove_if (dist.begin(), dist.end(), ptr_elim),
// 		 dist.end());    
//   }
  
  t_elim.stop();
}

void CandidateTeeth::get_range(const int root, const ToothBody &s,
			       IntPair &range,
			       const vector<vector<int>> &zones)
{
  range = IntPair(-1, -1);
  int start = zones[root][s.start], end = zones[root][s.end];

  if (start == end) {//if same endpoints
    if (start < 0) {//singleton
      range = IntPair(-start, -start);
    }//else empty
    return;
  }

  //now different endpoints
  if (start < 0) {
    range.first = -start; //first definitely start
    range.second = fabs(end);
    return;
  }

  //now diff endpoints, start >= 0
  if (end < 0) {
    range.second = -end; //second definitely end
    if (start == -end)
      range.first = start; //start is itself
    else
      range.first = start + 1;// fmin(start + 1, -end) start is at most end
    return;
  }

  //now diff endpoints, both >= 0
  range = IntPair(start + 1, end);// fmax(start + 1, end - 1));
}

bool CandidateTeeth::root_equivalent(const int root, const ToothBody &s1,
				     const ToothBody &s2,
				     const vector<vector<int>> &zones)
{
  IntPair s1_range, s2_range;
    
  get_range(root, s1, s1_range, zones);
  get_range(root, s2, s2_range, zones);
  return s1_range == s2_range;
}

bool CandidateTeeth::root_equivalent(const int root, const ToothBody &s1,
				     const ToothBody &s2) const
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
  if ((cut_end - cut_start + 1) > (arg->G_s.node_count - 2))
      return 0;
  
  //distant declarations
  vector<int>  &marks = arg->node_marks;
  vector<int> &tour = arg->tour_nodes;
  vector<int> &perm = arg->perm;
  
  SupportGraph &G = arg->G_s;
  
  std::unordered_map<int, double> &rb_sums = arg->rb_sums;
  
  int ncount = G.node_count;
  
  double rb_lower = cut_val - (1.5 - Epsilon::Cut);

  //right-adjacent declarations
  ToothBody &old_seg = arg->old_seg;
  vector<vector<int>> &zones = arg->adj_zones;
  vector<util::SquareUT<ToothList::reverse_iterator>> &ranges = arg->ranges;

  //distant add
  if (cut_start == old_seg.start) {
    for (int i = old_seg.end + 1; i <= cut_end; ++i) {
      marks[i] = 1;
      rb_sums.erase(i);
    }

    for (int i = old_seg.end + 1; i <= cut_end; ++i) {
      SNode vx = G.nodelist[tour[i]];
      for (int k = 0; k < vx.s_degree; ++k) {
	int root_perm = perm[vx.adj_objs[k].other_end];

	if (marks[root_perm] == 0)
	  rb_sums[root_perm] += vx.adj_objs[k].lp_weight;
      }
    }
  } else {//on a new degree node
    //clean up after the old one
    for (int i = 0; i < ncount; ++i) {
      marks[i] = 0;
    }
    rb_sums.clear();

    //set up for the new one
    marks[cut_start] = 1;

    SNode vx = G.nodelist[tour[cut_start]];
    for (int k = 0; k < vx.s_degree; ++k) {
      int root_perm = perm[vx.adj_objs[k].other_end];

      if (marks[root_perm] == 0)
	rb_sums[root_perm] += vx.adj_objs[k].lp_weight;
    }    
  }

  for (auto &kv : rb_sums) {
    int root = kv.first;
    double rb_sum = kv.second;
    if ((cut_start == cut_end) && (root > cut_end) && (root != (ncount - 1)))
      continue;
    if (((cut_end - cut_start + 1) == (ncount - 2)) && (root > cut_end))
      continue;

    if (rb_sum > rb_lower) {
        ToothList &teeth = arg->light_teeth[root];
        double abs_slack = fabs(cut_val - rb_sum - 1);
        double new_slack = (abs_slack < Epsilon::Zero) ? 0 : abs_slack;

      try {
          add_tooth(teeth, zones, ranges, root, cut_start, cut_end, new_slack);
      } catch (...) { CMR_SET_GOTO(rval, "Couldn't push back dist tooth. "); }
    }
  }


 CLEANUP:
  if (rval)
    cerr << "CandidateTeeth::tooth_cb failed.\n";
  old_seg = ToothBody(cut_start, cut_end, slack); //right adjacent update
  return rval;
}

inline void CandidateTeeth::add_tooth(ToothList &teeth,
                                      const vector<vector<int>> &zones,
                                      vector<
                                      util::SquareUT<ToothList::reverse_iterator>>
                                      & ranges,
                                      const int root, const int body_start,
                                      const int body_end, const double slack)
{
    bool elim = false;
    ToothBody body(body_start, body_end);
    IntPair range;
    get_range(root, body, range, zones);
  
    if (!teeth.empty()) {        
        ToothList::reverse_iterator &rit = ranges[root](range.first,
                                                        range.second);

        if (rit != teeth.rend()) {
            elim = true;
            double old_slack = (*rit)->slack;
            if (slack < old_slack)
                *rit = util::make_unique<SimpleTooth>(root, body, slack);
        }
    }

  if (!elim) {
      teeth.emplace_back(util::make_unique<SimpleTooth>(root, body, slack));
      ranges[root](range.first, range.second) = teeth.rbegin();
  }
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

  if (T.body_start == T.body_end) {
    cout << "}) -- slack " << T.slack << "\n";
    return;
  }

  if (!full) {
    cout << ", ..., " << bt[T.body_end] << "}) -- slack " << T.slack << "\n";
    return;
  }

  bool comma_sep = T.body_end - T.body_start <= 20;
  int i = T.body_start;
  
  cout << ", ";
  while (i++ != T.body_end) {
    cout << bt[i];
    if (i == T.body_end)
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
  t_comp_elim.report(true);
  t_sort.report(true);
  t_all.report(true);
}

}
