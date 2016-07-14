#include<iomanip>

#include "tooth.h"
#include "PSEP_util.h"

using namespace std;

int PSEP_CandTooth::SimpleTooth::ncount;
SupportGraph * PSEP_CandTooth::SimpleTooth::G_s;
int * PSEP_CandTooth::SimpleTooth::edge_marks;
int * PSEP_CandTooth::SimpleTooth::best_tour_nodes;

void PSEP_CandTooth::build_collection(){
  SimpleTooth::edge_marks = &edge_marks[0];

  for(int i = 0; i < SimpleTooth::ncount; i++){
    if(!light_teeth[i].empty())
      light_teeth[i].clear();
    if(!heavy_teeth[i].empty())
      heavy_teeth[i].clear();
    
    find_root_adjacent_teeth(i);
    find_root_distant_teeth(i);

    
    if(!light_teeth[i].empty())
      light_teeth[i].sort(SimpleTooth::p_greater);
    if(!heavy_teeth[i].empty())
      heavy_teeth[i].sort(SimpleTooth::p_greater);
  }
}

void PSEP_CandTooth::find_root_adjacent_teeth(const int root){
  int ncount = SimpleTooth::ncount;
  int body_start = (root + 1) % ncount, body_end;
  double lhs = 0.0;
  int rhs = -1;

  bool found_dup;

  for(int i = 1; i < ncount - 1; i++){
    found_dup = false;
    body_end = (root + i) % ncount;
    std::unique_ptr<SimpleTooth>
      cand(new SimpleTooth(root, body_start, body_end));
    int new_vx = best_tour_nodes[body_end];

    cand->increment_slack(new_vx, &lhs, &rhs);

    if(cand->slack >= 1 - LP::EPSILON || cand->slack < 0)
      continue;

    if(cand->body_size > (ncount - 2) / 2)
      cand->complement();

    if(cand->slack < 0.5){
      if(cand->body_size == 1){
	if(!light_teeth[cand->body_start].empty())
	  for(list<std::unique_ptr<SimpleTooth> >::reverse_iterator
		orig = light_teeth[cand->body_start].rbegin();
	      orig != light_teeth[cand->body_start].rend(); orig++){
	    if((*orig)->body_size > 1) break;
	    if((*orig)->body_start == cand->root &&
	       cand->body_start == (*orig)->root){
	      found_dup = true;
	      break;
	    }
	  }
	if(found_dup) continue;
      }      
      light_teeth[root].push_back(std::move(cand));
    } else{
      if(cand->body_size == 1){
	if(!heavy_teeth[cand->body_start].empty())
	  for(list<std::unique_ptr<SimpleTooth> >::reverse_iterator
		orig = heavy_teeth[cand->body_start].rbegin();
	      orig != heavy_teeth[cand->body_start].rend(); orig++){
	    if((*orig)->body_size > 1) break;
	    if((*orig)->body_start == cand->root &&
	       cand->body_start == (*orig)->root){
	      found_dup = true;
	      break;
	    }
	  }
	if(found_dup) continue;
      }
      heavy_teeth[root].push_back(std::move(cand));
    }
  }

  for(int k = 0; k < ncount; k++) edge_marks[k] = 0;
}


void PSEP_CandTooth::find_root_distant_teeth(const int root){
  int ncount = SimpleTooth::ncount;
  int body_start, body_end;
  double lhs;
  int rhs;

  bool found_dup;

  for(int i = 2; i < ncount - 1; i++){
    body_start = (root + i) % ncount;
    for(int k = 0; k < ncount; k++) edge_marks[k] = 0; //REPLACE???
    lhs = 0.0; rhs = -1;
    for(int j = i; j < ncount - 1; j++){
      body_end = (root + j) % ncount;
      unique_ptr<SimpleTooth> cand(new SimpleTooth(root, body_start, body_end));
      int new_vx = best_tour_nodes[body_end];

      cand->increment_slack(new_vx, &lhs, &rhs);

      if(cand->slack >= 1 - LP::EPSILON || cand->slack < 0 )
	continue;

      if(cand->body_size > (ncount - 2) / 2)
	cand->complement();

      if(cand->slack < 0.5){
	if(cand->body_size == 1){
	  if(!light_teeth[cand->body_start].empty())
	    for(list<unique_ptr<SimpleTooth> >::reverse_iterator
		  orig = light_teeth[cand->body_start].rbegin();
		orig != light_teeth[cand->body_start].rend(); orig++){
	      if((*orig)->body_size > 1) break;
	      if((*orig)->body_start == cand->root &&
		 cand->body_start == (*orig)->root){
		found_dup = true;
		break;
	      }
	    }
	  if(found_dup) continue;
	}      
	light_teeth[root].push_back(std::move(cand));
      } else {
	if(cand->body_size == 1){
	  if(!heavy_teeth[cand->body_start].empty())
	    for(list<unique_ptr<SimpleTooth> >::reverse_iterator
		  orig = heavy_teeth[cand->body_start].rbegin();
		orig != heavy_teeth[cand->body_start].rend(); orig++){
	      if((*orig)->body_size > 1) break;
	      if((*orig)->body_start == cand->root &&
		 cand->body_start == (*orig)->root){
		found_dup = true;
		break;
	      }
	    }
	  if(found_dup) continue;
	}
	heavy_teeth[root].push_back(std::move(cand));
      }
    }
  }

  for(int k = 0; k < ncount; k++) edge_marks[k] = 0;
}


bool PSEP_CandTooth::SimpleTooth::body_contains(const int perm_node){
  if(perm_node == root)
    return false;

  return (body_start <= body_end) ?
    (body_start <= perm_node && perm_node <= body_end) :
    (body_start <= perm_node || perm_node <= body_end);
}

//returns whether bodies T C R
bool PSEP_CandTooth::SimpleTooth::C_body_subset(const SimpleTooth &T,
						const SimpleTooth &R){
  if(T.root == R.root){
    if(R.body_start <= R.body_end){ //if R is contiguous
      if(T.body_start > T.body_end){//R cannot contain a sandwich tooth T
	return false;
      } else {
	return (R.body_start <= T.body_start && T.body_end <= R.body_end);
      }
    } else { //if R wraps
      return ( (R.body_start <= T.body_start &&
		(R.body_start <= T.body_end || T.body_end <= R.body_end)) ||
	       (T.body_start <= R.body_end && T.body_end <= R.body_end));
    }
  } else {
    std::cerr << "TRIED TO USE TEETH W DIF ROOT!!!\n";
    exit(1);
  }
}

void PSEP_CandTooth::SimpleTooth::complement(){
  int c_start, c_end;
  bool c_sandwich;

  if((body_start == ((root + 1) % ncount)) ||
     (root == ((body_end + 1) % ncount))){
    c_sandwich = false;
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
    c_sandwich = !(sandwich);
  }

  body_start = c_start;
  body_end = c_end;
  sandwich = c_sandwich;
}

//new vx is some actual best_tour_nodes[k], not k
void PSEP_CandTooth::SimpleTooth::increment_slack(const int new_vx,
						  double *lhs_p, int *rhs_p){
  SNode *current_node = &(G_s->nodelist[new_vx]);
  int other_end;
  double lp_weight;

  *rhs_p += 2;
  edge_marks[new_vx] = 1;

  for(int i = 0; i < current_node->s_degree; i++){
    other_end = current_node->adj_objs[i].other_end;
    lp_weight = current_node->adj_objs[i].lp_weight;

    if(edge_marks[other_end] == 1){
      *lhs_p += (2 * lp_weight);
      continue;
    }

    if(other_end == root)
      *lhs_p += lp_weight;
  }

  slack = *rhs_p - *lhs_p;
}

void PSEP_CandTooth::SimpleTooth::print(){
  int ncount = SimpleTooth::ncount;
  int current_node = body_start;
  int upper_limit = body_size + ((int) sandwich);

  cout << "Root: " << best_tour_nodes[root];
  cout << ", Body size: " << body_size << "\n";
  cout << "Body: " << best_tour_nodes[body_start] << "\n";
  for(int i = 1; i < upper_limit; i++){
    current_node = best_tour_nodes[(body_start + i) % ncount];
    if(current_node != best_tour_nodes[root])
      cout << setw(7) << current_node << "\n";
  }
  cout << "Slack: " << slack << "\n";
  cout << "Sandwich: " << sandwich << "\n";
}
