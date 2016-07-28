#include<iostream>
#include<math.h>

#include "LPfixing.h"

using namespace std;

int PSEP_LPfix::price(){
  int rval = 0;
  int fixtotal = 0;
  int ecount = m_graph.edge_count;
  vector<double> redcosts(ecount);
  double GAP, lp_LB, cur_red, cur_abs;
  double start;

  int infeasible = 0;
  start = PSEP_zeit();
  rval = PSEPlp_primal_opt(&m_lp, &infeasible);
  if(rval) goto CLEANUP;
  start = PSEP_zeit() - start;
  cout << "Optimized LP relaxation in "
       << start << "s, objective val: ";

  rval = PSEPlp_objval(&m_lp, &lp_LB);
  if(rval) goto CLEANUP;
  cout << lp_LB << "\n";

  GAP = m_min_tour_value - lp_LB;
  cout << "Integrality gap: " << GAP << "\n";

  rval = PSEPlp_get_redcosts(&m_lp, &redcosts[0]);
  if(rval) goto CLEANUP;
  
  for(int i = 0; i < ecount; i++){
    cur_red = redcosts[i];
    cur_abs = fabs(cur_red);
    if(cur_abs <= GAP - 1){
      edge_delset[i] = FixStats::LEAVE;
      continue;
    }

    fixtotal++;
    cout << "Edge " << i << "has reduced cost " << redcosts[i] << ", ";
    if(cur_red > 0){
      edge_delset[i] = FixStats::DELETE;
      cout << "deleting it.\n";
      continue;
    }

    if(cur_red < 0){
      edge_delset[i] = FixStats::FIXED;
      cout << "permanently fixing to 1.\n";
    }
  }

  cout << fixtotal << " total edges can be fixed\n";
  

 CLEANUP:
  if(rval)
    cerr << "Error entry point: LPfixing::price() \n";
  return rval;
}
