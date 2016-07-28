#include<iostream>
#include<iomanip>
#include<algorithm>
#include<math.h>

#include "LPfixing.h"

using namespace std;

int PSEP_LPfix::price(int *clamptotal, int *deltotal){
  int rval = 0;
  int ecount = m_graph.edge_count;
  vector<double> redcosts(ecount);
  double GAP, lp_LB, cur_red, cur_abs;
  double opt_time;
  double total_time = PSEP_zeit();
  *clamptotal = 0; *deltotal = 0;

  edge_delset.resize(ecount);

  int infeasible = 0;
  opt_time = PSEP_zeit();
  rval = PSEPlp_primal_opt(&m_lp, &infeasible);
  if(rval) goto CLEANUP;
  opt_time = PSEP_zeit() - opt_time;
  cout << "Optimized LP relaxation in "
       << setprecision(0) << opt_time << "s, obj val: " << setprecision(6);

  rval = PSEPlp_objval(&m_lp, &lp_LB);
  if(rval) goto CLEANUP;
  cout << lp_LB << ", ";

  GAP = m_min_tour_value - lp_LB;
  cout << "Integrality gap: " << GAP << ".\n";

  rval = PSEPlp_get_redcosts(&m_lp, &redcosts[0]);
  if(rval) goto CLEANUP;
  
  for(int i = 0; i < ecount; i++){
    cur_red = redcosts[i];
    cur_abs = fabs(cur_red);
    if(cur_abs <= GAP - 1){
      edge_delset[i] = FixStats::LEAVE;
      continue;
    }

    if(cur_red > 0){
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
    cout << *clamptotal << " edges shall be fixed to 1, ";
  if(*deltotal != 0)
    cout << *deltotal << " edges shall be deleted, "
	 << "computed in " << setprecision(0)
	 << (PSEP_zeit() - total_time) << "s\n";
  

 CLEANUP:
  if(rval)
    cerr << "Error entry point: LPfixing::price() \n";
  return rval;
}

int PSEP_LPfix::fixup(){
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

int PSEP_LPfix::delete_cols(){
  int rval = PSEPlp_delsetcols(&m_lp, &edge_delset[0]);
  if(rval)
    cerr << "Error entry point: LPfix::delete_cols \n";
  return rval;
}

void PSEP_LPfix::delete_edges(){
  int orig_ecount = m_graph.edge_count;
  double total_time = PSEP_zeit();
  for(int i = 0; i < m_graph.edge_count; i++)
    edge_lookup[IntPair(edges[i].end[0], edges[i].end[1])] = edge_delset[i];

  vector<Edge> new_graph_edges;
  vector<int> new_best_edges;

  for(int i = 0; i < edge_delset.size(); i++){
    if(edge_delset[i] != -1){
      new_best_edges.push_back(best_tour_edges[i]);
      new_graph_edges.emplace_back(edges[i]);
    }
  }

  edges.clear();
  edges = new_graph_edges;
  best_tour_edges.clear();
  best_tour_edges = new_best_edges;
  
  m_graph.edge_count = edges.size();
  delta.resize(m_graph.edge_count);
  m_lp_edges.resize(m_graph.edge_count);

  total_time = PSEP_zeit() - total_time;

  cout << "Number of edges is now " << edges.size() << ", deleted in "
       << total_time << "s.\n (Down from " << orig_ecount << ", ";
  cout << "ratio "
       << setprecision(3) << ((double) edges.size() / m_graph.node_count)
       << " * ncount)\n";
}


int PSEP_LPfix::redcost_fixing(){
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
  return rval;
}
