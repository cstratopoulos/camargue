#include<cmath>

#include "BButils.h"

using namespace std;
using namespace PSEP::BB;

int Constraints::enforce(unique_ptr<TreeNode> &v){
  if(v->type() == TreeNode::NType::ROOT){
    cerr << "BB::Constraints::enforce tried to enforce ROOT\n";
    return 1;
  }

  if(v->status() != TreeNode::NStat::UNVISITED){
    cerr << "BB::Constraints::enforce tried to enforce visited node\n";
    return 1;
  }

  int clamp_edge = v->branch_edge;
  int rval = 0;

  if(v->type() == TreeNode::NType::LEFT){
    rval = add_left_clamp(clamp_edge);
    if(rval) goto CLEANUP;
  } else {
    rval = add_right_branch(clamp_edge);
    if(rval) goto CLEANUP;
  }

 CLEANUP:
  if(rval)
    cerr << "Entry point: Constraints::enforce\n";
  return rval;  
}

int Constraints::compute_branch_edge(){
  int branch_edge = -1, index;
  double max_under_half = 0.0, min_over_half = 1.0;
  double LB = 0.0, UB = 1.0;
  double lp_weight;
  int max_len = 0;

  for(int i = 0; i < support_indices.size(); i++){
    index = support_indices[i];
    lp_weight = m_lp_edges[index];
    if(fabs(lp_weight - 0.5) < LP::EPSILON){
      max_under_half = 0.5; min_over_half = 0.5;
      break;
    }

    if(lp_weight < 0.5){
      if(lp_weight > max_under_half)
	max_under_half = lp_weight;
      continue;
    }

    if(lp_weight > 0.5){
      if(lp_weight < min_over_half)
	min_over_half = lp_weight;
      continue;
    }
  }

  LB = 0.75 * max_under_half;
  UB = min_over_half + (0.25 * (1.0 - min_over_half));

  for(int i = 0; i < support_indices.size(); i++){
    index = support_indices[i];
    lp_weight = m_lp_edges[index];
    if(lp_weight < LB || lp_weight > UB)
      continue;
    
    if(edges[index].len > max_len){
      branch_edge = index;
      max_len = edges[index].len;
    }
  }

  return branch_edge;
}

int Constraints::add_left_clamp(const int edge){
  double bound = best_tour_edges[edge];
  char lower_or_upper = (bound == 0.0) ? 'U' : 'L';

  int rval = PSEPlp_clampbnd(&m_lp, edge, lower_or_upper, bound);
  if(rval)
    cerr << "Problem in BB::Constraints::add_left_clamp\n";

  return rval;
}

int Constraints::remove_left_clamp(const int edge){
  double bound = !best_tour_edges[edge];
  char lower_or_upper = (bound == 0.0) ? 'L' : 'U';

  int rval = PSEPlp_clampbnd(&m_lp, edge, lower_or_upper, bound);
  if(rval)
    cerr << "Problem in BB::Constraints::remove_left_clamp\n";

  return rval;
}

inline int Constraints::add_right_branch(const int edge){
  return RBranch.active() ? explore_right(edge) :
    add_first_right_rows(edge);
};

int Constraints::explore_right(const int edge){
  int rval = 0, rownum;
  std::map<int, int>::iterator row_it =
    RBranch.edge_row_lookup.find(edge);

  if(row_it == RBranch.edge_row_lookup.end()){
    cerr << "Tried to explore right branch for invalid constraint\n"; return 1;
  }

  rownum = row_it->second;

  char sense = 'E';

  rval = PSEPlp_chgsense(&m_lp, 1, &edge, &sense);
  if(rval)
    cerr << "Problem in BB::Constraints::explore_right\n";

  return rval;
}

int Constraints::add_first_right_rows(const int edge){
  if(RBranch.active()){
    cerr << "Tried to add new right rows but right branch was already active\n";
    return 1;
  }

  int range_start = PSEPlp_numrows(&m_lp), range_end;
  int current_rownum = range_start;
  int newnz = 2, newrows = 1, rmatbeg = 0;
  std::array<int, 2> rmatind = {edge, -1};
  std::array<double, 2> rmatval;
  double RHS;
  char sense = 'L';

  for(int partner = 0; partner < PSEPlp_numcols(&m_lp); partner++){
    if(EdgeStats.Left.count(partner) != 0 ||
       EdgeStats.Right.count(partner) != 0 || partner == edge)
      continue;

    rmatind[1] = partner;
    compute_right_row(edge, partner, rmatval, RHS);

    if(PSEPlp_addrows(&m_lp, newrows, newnz, &RHS, &sense,
		      &rmatbeg, &rmatind[0], &rmatval[0])){
      cerr << "Problem in BB::Constraints::add_first_right_rows\n";
      return 1;
    }

    RBranch.add_edge_row(partner, current_rownum++);
  }

  RBranch.first_right = edge;
  range_end = PSEPlp_numrows(&m_lp) - 1;
  RBranch.constraint_range = IntPair(range_start, range_end);

  return 0;
}


