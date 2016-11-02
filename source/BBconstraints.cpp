#include<algorithm>

#include<cmath>

#include "BBconstraints.hpp"

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
    cout << "\n  Added left clamp on edge " << clamp_edge << "\n";
  } else {
    rval = add_right_branch(clamp_edge);
    if(rval) goto CLEANUP;
    cout << "\n  Added right branch on edge " << clamp_edge << "\n";
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
    cout << "Removed left clamp on edge " << clamp_edge << "\n";
  } else {
    rval = remove_right(clamp_edge);
    if(rval) goto CLEANUP;
    cout << "Removed right branch constraint(s) on edge "
	 << clamp_edge;
    cout << ". first right: " << RBranch.first_right_edge() << ", "
	 << "now active: " << RBranch.active() << "\n";
  }

  EdgeStats.remove_branch_var(v);

 CLEANUP:
  if(rval)
    cerr << "Problem in BB::Constraints::unenforce\n";
  return rval;
}

bool Constraints::naive_compatible(const int clamp, const int partner){
  int best_part  = best_tour_edges[partner];
  double lp_clamp = m_lp_edges[clamp], lp_part = m_lp_edges[partner];

  switch(best_part){
  case 0:
    return lp_clamp < lp_part;
  case 1:
  default:
    return lp_part < lp_clamp;
  }
}

int Constraints::compute_branch_edge(){
  int branch_edge = -1, index;
  double max_under_half = 0.0, min_over_half = 1.0;
  double LB = 0.0, UB = 1.0;
  double lp_weight;
  int max_len = 0;
  naive_branch_candidates.clear();
  
  for(int i = 0; i < support_indices.size(); i++){
    index = support_indices[i];
    if(EdgeStats.Right.count(index) || EdgeStats.Left.count(index)) continue;
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
    if(lp_weight < LB || lp_weight > UB || EdgeStats.Right.count(index)
       || EdgeStats.Left.count(index))
      continue;

    if(Strategy == BranchPlan::Naive)
      naive_branch_candidates.push_back(index);
    else
      if(edges[index].len > max_len){
	branch_edge = index;
	max_len = edges[index].len;
      }
  }

  if(Strategy == BranchPlan::Naive){
    sort(naive_branch_candidates.begin(), naive_branch_candidates.end(),
	 [this](const int a, const int b) -> bool {
	   return edges[a].len > edges[b].len; });
    for(int i = 0; i < naive_branch_candidates.size(); i++){
      int edge = naive_branch_candidates[i];
      for(int partner = 0; partner < m_lp_edges.size(); partner++){
	if(naive_compatible(edge, partner)){
	  naive_edge_partner = partner;
	  return branch_edge;
	}
      }
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

void Constraints::compute_right_row(const int clamp, const int partner,
				   std::array<double, 2> &rmatval, double &RHS){
  double clamp_best = best_tour_edges[clamp],
    partner_best = best_tour_edges[partner];
  RHS = clamp_best - partner_best;
  rmatval[0] = 2.0 * clamp_best - 1;
  rmatval[1] = 1 - 2.0 * partner_best;
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
  rmatval[0] = 2 * clamp_best - 1;
  rmatval[1] = 1 - 2 * partner_best;

  return 0;
}

int Constraints::add_right_branch(const int edge){
  int rval = 0;

  switch(Strategy){
  case BranchPlan::Main:
    rval = RBranch.active() ? explore_right(edge) :
      add_main_right_rows(edge);
    break;
  case BranchPlan::Naive:
    rval = 1;
    cout << "Alternate case goes here\n";
    break;
  }

  if(rval)
    cerr << "Problem in Constraints::add_right_branch\n";
  return rval;
}

int Constraints::explore_right(const int edge){
  int rval = 0, rownum;
  std::map<int, int>::iterator row_it =
    RBranch.edge_row_lookup.find(edge);

  if(row_it == RBranch.edge_row_lookup.end()){
    cerr << "Tried to explore right branch for invalid constraint "
	 << "on edge " << edge << "\n";
    return 1;
  }

  rownum = row_it->second;

  char sense = 'E';

  rval = PSEPlp_chgsense(&m_lp, 1, &rownum, &sense);
  if(rval)
    cerr << "Problem in BB::Constraints::explore_right\n";

  return rval;
}

int Constraints::add_main_right_rows(const int edge){
  if(RBranch.active()){
    cerr << "Problem in Constraints::add_main_right_rows: right branch "
	 << "was already active\n";
    return 1;
  }

  if(LPCore.rebuild_basis(true)){
    cerr << "Problem in Constraints::add_main_right_rows\n";
  }

  int range_start = PSEPlp_numrows(&m_lp), range_end;
  int numrows = range_start; 
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
      cerr << "Problem in BB::Constraints::add_main_right_rows\n";
      return 1;
    }

    RBranch.add_edge_row(partner, current_rownum++);
  }

  RBranch.first_right = edge;
  range_end = PSEPlp_numrows(&m_lp) - 1;
  RBranch.constraint_range = IntPair(range_start, range_end);

  cout << "  Added " << (PSEPlp_numrows(&m_lp) - numrows) << " rows for "
       << " the right branch";
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
  }

  RBranch.first_right = -1;
  RBranch.constraint_range = IntPair(0,0);

  return 0;
}

int Constraints::prune(){
  int rval = 0;
  vector<int> delset;
  IntPair skiprange = RBranch.skiprange;
  int num_removed;

  if(RBranch.active()){
    rval = LPCore.rebuild_basis(num_removed, skiprange, delset);
    if(rval) goto CLEANUP;

    if(num_removed > 0){
      rval = RBranch.update_range(delset);
      if(rval) goto CLEANUP;
    }
  } else {
    rval = LPCore.rebuild_basis(true);
    if(rval) goto CLEANUP;
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in Constraints::prune\n";
  return rval;
}


