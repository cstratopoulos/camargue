#include "tooth.hpp"

#include <iostream>
#include <algorithm>

using namespace PSEP;
using namespace std;

CandidateTeeth::CandidateTeeth(vector<int> &_edge_marks,
			       vector<int> &_best_tour_nodes,
			       SupportGraph &_G_s) :
    edge_marks(_edge_marks), best_tour_nodes(_best_tour_nodes), G_s(_G_s)
{
  light_teeth.resize(best_tour_nodes.size());
}

int CandidateTeeth::get_light_teeth()
{
  int rval = 0;
  clear_collection();

  double ft = 0, ft_d = 0, st = 0;

  for(int i = 0; i < G_s.node_count; i++){
    //    cout << "=== ROOT " << i << ", ";
      
    double ft_i = zeit();
    rval = get_adjacent_teeth(i);
    if(rval) goto CLEANUP;
    ft += zeit() - ft_i;

    //    cout << light_teeth[i].size() << " adj teeth, ";

    double ft_d_i = zeit();
    rval = get_distant_teeth(i);
    if(rval) goto CLEANUP;
    ft_d += zeit() - ft_d_i;

    //    cout << light_teeth[i].size() << " total after dist\n";// << endl;

    double st_i = zeit();
    if(!light_teeth[i].empty()) //lambda to sort by decreasing body size
      sort(light_teeth[i].begin(), light_teeth[i].end(),
	   [this](const SimpleTooth::Ptr &T,
		  const SimpleTooth::Ptr &R) -> bool {
	     return body_size(*T) > body_size(*R);
	   });
    st += zeit() - st_i;
  }

  cout << ft << "s finding adjacent teeth, " << ft_d
       << "s finding distant teeth, " << st << "s sorting" << endl;

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
    SimpleTooth::Ptr cand(new SimpleTooth(root, body_start, body_end));

    increment_slack(*cand, new_vx, lhs, rhs);

    //NOTE: we are completely ignoring heavy teeth for now
    if(cand->slack >= 0.5 - LP::EPSILON || cand-> slack < 0)
      continue;

    if(body_size(*cand) > (ncount - 2) / 2)
      complement(*cand);

    if(body_size(*cand) == 1){
      if(cand->root > cand->body_start &&
	 !light_teeth[cand->body_start].empty()){
	bool found_dup = false;
	for(auto orig = light_teeth[cand->body_start].rbegin();
	    orig != light_teeth[cand->body_start].rend(); orig++){
	  if(body_size(**orig) > 1) break;

	  if((*orig)->body_start == cand->root &&
	     cand->body_start == (*orig)->root){
	    found_dup = true;
	    break;
	  }
	}
	if(found_dup)
	  continue;
      }
    }

    try { light_teeth[root].push_back(move(cand)); } catch(...){
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
