#include<cmath>

#include "BBconstraints.h"

using namespace std;
using namespace PSEP::BB;

int Constraints::enforce(unique_ptr<TreeNode> &v){
  if(v->type() == NodeType::ROOT){
    cout << "BB::Constraints::enforce is considering ROOT\n";
    return 0;
  }

  if(v->status() != NodeStat::UNVISITED){
    cerr << "BB::Constraints::enforce tried to enforce visited node\n";
    return 1;
  }

  int clamp_edge = v->edge();
  int rval = 0;

  if(v->type() == NodeType::LEFT){
    rval = add_left_clamp(clamp_edge);
    if(rval) goto CLEANUP;
  } else {
    rval = add_right_branch(clamp_edge);
    if(rval) goto CLEANUP;
  }

  EdgeStats.add_branch_var(v);

 CLEANUP:
  if(rval)
    cerr << "Entry point: BB::Constraints::enforce\n";
  return rval;  
}

int Constraints::unenforce(unique_ptr<TreeNode> &v){
  if(v->type() == NodeType::ROOT){
    cout << "BB::Constraints::unenforce is consideering ROOT\n";
    return 0;
  }

  if(v->status() == NodeStat::UNVISITED){
    cerr << "BB::Constraints::unenforce called on an univisited node\n";
    return 1;
  }

  int clamp_edge = v->edge();
  int rval = 0;

  if(v->type() == NodeType::LEFT){
    rval = remove_left_clamp(clamp_edge);
    if(rval) goto CLEANUP;
  } else {
    rval = remove_right(clamp_edge);
    if(rval) goto CLEANUP;
  }

  EdgeStats.remove_branch_var(v);

 CLEANUP:
  if(rval)
    cerr << "Problem in BB::Constraints::unenforce\n";
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
  double bound = 1 - best_tour_edges[edge];
  char lower_or_upper = (bound == 0.0) ? 'L' : 'U';

  int rval = PSEPlp_clampbnd(&m_lp, edge, lower_or_upper, bound);
  if(rval)
    cerr << "Problem in BB::Constraints::remove_left_clamp\n";

  return rval;
}

int Constraints::compute_right_update(const int clamp, const int partner,
				      array<double, 2> &rmatval, double &RHS,
				      const vector<double> &new_tour){
  double clamp_best, partner_best;

  if(fabs(new_tour[clamp]) < LP::EPSILON)
    clamp_best = 0;
  else if (fabs(new_tour[clamp]) > 1 - LP::EPSILON)
    clamp_best = 1;
  else {
    cerr << "New tour is not integral\n"; return 1;
  }

  if(fabs(new_tour[partner]) < LP::EPSILON)
    partner_best = 0;
  else if (fabs(new_tour[partner]) > 1 - LP::EPSILON)
    partner_best = 1;
  else {
    cerr << "New tour is not integral\n"; return 1;
  }

  RHS = clamp_best - partner_best;
  rmatval = {2 * clamp_best - 1, 1 - 2 * partner_best};

  return 0;
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

  rval = PSEPlp_chgsense(&m_lp, 1, &rownum, &sense);
  if(rval)
    cerr << "Problem in BB::Constraints::explore_right\n";

  return rval;
}

int Constraints::add_first_right_rows(const int edge){
  if(RBranch.active()){
    cerr << "Tried to add new right rows but right branch "
	 << "was already active\n";
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
       EdgeStats.FixedUp.count(partner) != 0 || partner == edge)
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

int Constraints::update_right_rows(){
  if(!RBranch.active())
    return 0;
  int rval = 0;
  int rownum, clamp = RBranch.first_right_edge(), partner;
  array<double, 2> rmatval;
  double RHS;
  vector<double> newtour(best_tour_edges.size());
  bool clamp_dif;

  rval = PSEPlp_x(&m_lp, &newtour[0]);
  if(rval) goto CLEANUP;

  clamp_dif = (fabs(best_tour_edges[clamp] - newtour[clamp]) >= LP::EPSILON);

  for(map<int, int>::const_iterator it = RBranch.edge_row_lookup.begin();
      it != RBranch.edge_row_lookup.end(); it++){
    partner = it->first;
    if(clamp_dif ||
       (fabs(best_tour_edges[partner] - newtour[partner]) >= LP::EPSILON)){
      rownum = it->second;
      rval = compute_right_update(clamp, partner, rmatval, RHS, newtour);
      if(rval) goto CLEANUP;

      rval = (PSEPlp_chgcoef(&m_lp, rownum, clamp, rmatval[0]) ||
	      PSEPlp_chgcoef(&m_lp, rownum, partner, rmatval[1]) ||
	      PSEPlp_chgcoef(&m_lp, rownum, -1, RHS));
      if(rval) goto CLEANUP;
    }
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in Constraints::update_right_rows\n";
  return rval;
}

int Constraints::remove_right(const int edge){
  if(!RBranch.active()){
    cerr << "Called Constraints::remove_right when none were active\n";
    return 1;
  }

  if(edge != RBranch.first_right_edge()){
    map<int, int>::const_iterator it = RBranch.edge_row_lookup.find(edge);
    if(it == RBranch.edge_row_lookup.end()){
      cerr << "Constraints::remove_right tried to undo nonexistant row\n";
      return 1;
    }

    int rownum = it->second;
    char sense = 'L';
    if(PSEPlp_chgsense(&m_lp, 1, &rownum, &sense)){
      cerr << "Constraints::remove_right couldn't switch sense\n";
      return 1;
    }

    return 0;
  }

  if(PSEPlp_delrows(&m_lp, RBranch.constraint_range.first,
		    RBranch.constraint_range.second)){
    cerr << "Constraints::remove_right failed to delete bunch of rows\n";
    return 1;
    RBranch.first_right = -1;
  }

  return 0;
}


