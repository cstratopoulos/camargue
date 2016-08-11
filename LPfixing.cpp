#include<iostream>
#include<iomanip>
#include<algorithm>
#include<math.h>

#include "LPfixing.h"

using namespace std;
using namespace PSEP::LP;

int EdgeFix::price(int *clamptotal, int *deltotal){
  int rval = 0;
  int ecount = m_graph.edge_count;
  vector<double> redcosts(ecount);
  double GAP, lp_LB, cur_red, cur_abs;
  double opt_time;
  *clamptotal = 0; *deltotal = 0;

  edge_delset.resize(ecount);

  int infeasible = 0;
  opt_time = PSEP_zeit();
  rval = PSEPlp_primal_opt(&m_lp, &infeasible);
  if(rval) goto CLEANUP;
  opt_time = PSEP_zeit() - opt_time;
  cout << "Optimized LP in "
       << setprecision(0) << opt_time << "s, " << setprecision(6);

  rval = PSEPlp_objval(&m_lp, &lp_LB);
  if(rval) goto CLEANUP;

  GAP = m_min_tour_value - lp_LB;
  cout << "integrality gap: " << GAP << ".\n";

  rval = PSEPlp_get_redcosts(&m_lp, &redcosts[0]);
  if(rval) goto CLEANUP;
  
  for(int i = 0; i < ecount; i++){
    cur_red = redcosts[i];
    cur_abs = fabs(cur_red);
    if(cur_abs <= GAP - 1){
      edge_delset[i] = FixStats::LEAVE;
      continue;
    }

    if(cur_red > 0 && best_tour_edges[i] != 1){
      edge_delset[i] = FixStats::DELETE;
      (*deltotal)++;
      continue;
    }

    if(cur_red < 0){
      edge_delset[i] = FixStats::FIXED;
      (*clamptotal)++;
    }
  }

  if(*clamptotal != 0)
    cout << *clamptotal << " edges fixable, ";
  if(*deltotal != 0)
    cout << *deltotal << " edges removable. ";
  

 CLEANUP:
  if(rval)
    cerr << "Error entry point: LPfixing::price() \n";
  return rval;
}

int EdgeFix::fixup(){
  int rval = 0;
  int total_time = PSEP_zeit();

  for(int i = 0; i < edge_delset.size(); i++){
    if(edge_delset[i] == FixStats::FIXED){
      rval = PSEPlp_clampbnd(&m_lp, i, 'L', 1.0);
      if(rval) goto CLEANUP;
      edge_delset[i] = FixStats::LEAVE;
    }
  }

  total_time += PSEP_zeit() - total_time;
  
 CLEANUP:
  if(rval)
    cerr << "Error entry point: LPfix::fixup \n";
  return rval;
}

int EdgeFix::delete_cols(){
  int rval = PSEPlp_delsetcols(&m_lp, &edge_delset[0]);
  if(rval)
    cerr << "Error entry point: LPfix::delete_cols \n";
  return rval;
}

void EdgeFix::delete_edges(){
  double total_time = PSEP_zeit();
  for(int i = 0; i < m_graph.edge_count; i++)
    edge_lookup[IntPair(edges[i].end[0], edges[i].end[1])] = edge_delset[i];

  vector<Edge> new_graph_edges;
  vector<int> new_best_edges;
  vector<double> new_lp_edges;

  for(int i = 0; i < edge_delset.size(); i++){
    if(edge_delset[i] != -1){
      new_best_edges.push_back(best_tour_edges[i]);
      new_graph_edges.emplace_back(edges[i]);
      new_lp_edges.push_back(m_lp_edges[i]);
    }
  }

  edges.clear();
  edges = new_graph_edges;
  best_tour_edges.clear();
  best_tour_edges = new_best_edges;
  
  m_graph.edge_count = edges.size();
  delta.resize(m_graph.edge_count);
  m_lp_edges = new_lp_edges;

  total_time = PSEP_zeit() - total_time;

  cout << edges.size() << " remain, "
       << "ratio "
       << setprecision(3) << ((double) edges.size() / m_graph.node_count)
       << " * ncount)\n";
}


int EdgeFix::redcost_fixing(){
  int rval = 0;
  int clamptotal = 0, deltotal = 0;

  rval = price(&clamptotal, &deltotal);
  if(rval) goto CLEANUP;

  if(clamptotal > 0){
    rval = fixup();
    if(rval) goto CLEANUP;
  }

  if(deltotal > 0){
    rval = delete_cols();
    if(rval) goto CLEANUP;
    delete_edges();
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in LPfix::redcost_fixing.\n";
  cout << setprecision(6);
  return rval;
}
