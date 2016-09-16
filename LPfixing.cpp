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
  vector<double> redcosts;
  double GAP, lp_LB, cur_red, cur_abs;
  double opt_time;
  int infeasible = 0;

  *clamptotal = 0; *deltotal = 0;

  //
  cout << "In EdgeFix::price....printing best tour\n";
  PSEP::print_vec_nonzero(best_tour_edges);
  //
  
  try {
    redcosts.resize(ecount);
    edge_delset.resize(ecount);
  } catch(const bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for redcosts/delset, ");
  }

  opt_time = zeit();
  rval = PSEPlp_primal_opt(&m_lp, &infeasible);
  if(rval) goto CLEANUP;
  opt_time = zeit() - opt_time;
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
      //
      cout << "Deleting edge " << i << "\n";
      //
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
  int total_time = zeit();

  for(int i = 0; i < edge_delset.size(); i++){
    if(edge_delset[i] == FixStats::FIXED){
      rval = PSEPlp_clampbnd(&m_lp, i, 'L', 1.0);
      if(rval) goto CLEANUP;
      edge_delset[i] = FixStats::LEAVE;
    }
  }

  total_time += zeit() - total_time;
  
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
  double total_time = zeit();
  for(int i = 0; i < m_graph.edge_count; i++)
    edge_lookup[IntPair(edges[i].end[0], edges[i].end[1])] = edge_delset[i];

  for(int i = 0; i < edge_delset.size(); i++)
    if(edge_delset[i] == -1){
      edges[i].removable = true;
      best_tour_edges[i] = -1;
      m_lp_edges[i] = -1;
    }

  best_tour_edges.erase(remove(best_tour_edges.begin(),
				  best_tour_edges.end(), -1),
			best_tour_edges.end());
  m_lp_edges.erase(remove(m_lp_edges.begin(), m_lp_edges.end(), -1),
		   m_lp_edges.end());
  edges.erase(remove_if(edges.begin(), edges.end(),
			[](Edge &e){ return e.removable; }),
	      edges.end());
  

  
  m_graph.edge_count = edges.size();
  delta.resize(m_graph.edge_count);

  total_time = zeit() - total_time;

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
