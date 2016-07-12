#include "tooth.h"
#include "PSEP_util.h"

using namespace std;

int PSEP_CandTooth::SimpleTooth::ncount;
SupportGraph * PSEP_CandTooth::SimpleTooth::G_s;
int * PSEP_CandTooth::SimpleTooth::edge_marks;


void PSEP_CandTooth::find_root_adjacent_teeth(const int root){
  int ncount = SimpleTooth::ncount;
  int body_start = (root + 1) % ncount, body_end;
  double lhs = 0.0;
  int rhs = -1;

  bool found_dup;

  for(int i = 1; i < ncount - 1; i++){
    found_dup = false;
    body_end = (root + i) % ncount;
    unique_ptr<SimpleTooth> cand(new SimpleTooth(root, body_start, body_end));
    int new_vx = best_tour_nodes[body_end];

    cand->increment_slack(new_vx, &lhs, &rhs);

    if(cand->slack >= 1 - LP::EPSILON || cand->slack < 0)
      continue;

    if(cand->body_size() > (ncount - 2) / 2)
      cand->complement();

    if(cand->body_size() == 1){
      for(list<unique_ptr<SimpleTooth> >::iterator
	    orig = light_teeth[cand->body_start].begin();
	  orig != light_teeth[cand->body_start].end(); orig++){
	if((*orig)->body_size() > 1) continue;
	if((*orig)->body_start == cand->root &&
	   cand->body_start == (*orig)->root){
	  found_dup = true;
	  break;
	}
      }
      if(found_dup) continue;
    }

    if(cand->slack < 0.5)
      light_teeth[root].push_back(std::move(cand));
    else
      heavy_teeth[root].push_back(std::move(cand));
  }
}

int PSEP_CandTooth::SimpleTooth::body_size(){
  int add_one = (int) (!sandwich); //total is + 1 if tooth is not sandwiched

  return (body_start <= body_end) ? (body_end - body_start + add_one) :
    ((ncount - body_start) + body_end + add_one);
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
