#include<iostream>

#include "LPprune.h"

using namespace std;
using namespace PSEP;

int LPPrune::prune_cuts(int &num_removed){
  int rval = 0;
  num_removed = 0;
  int rowcount = PSEPlp_numrows(&m_lp);
  vector<int> delset(rowcount, 0);
  vector<double> slacks(rowcount - node_count, 0);

  rval = PSEPlp_getslack(&m_lp, &slacks[0], node_count, rowcount - 1);
  if(rval){
    cerr << "Problem in LPPrune::prune_cuts with getslack\n";
    return 1;
  }

  for(int i = 0; i < slacks.size(); i++){
    if(fabs(slacks[i]) >= LP::EPSILON){
      delset[node_count + i] = 1;
      num_removed++;
    }
  }

  rval = PSEPlp_delsetrows(&m_lp, &delset[0]);
  if(rval){
    cerr << "Problem in LPPrune::prune_cuts with delsetrows\n";
    return 1;
  }

  return 0;
}

int LPPrune::prune_with_skip(int &num_removed, IntPair skiprange,
			     vector<int> &delset){
  int rval = 0;
  num_removed = 0;
  int rowcount = PSEPlp_numrows(&m_lp);
  delset.clear();
  delset.resize(rowcount, 0);
  int skip_start = skiprange.first;
  int skip_end = skiprange.second;
  vector<double> slacks(rowcount - node_count, 0);

  rval = PSEPlp_getslack(&m_lp, &slacks[0], node_count, rowcount - 1);
  if(rval){
    cerr << "Problem in LPPrune::prune_with_skip with getslack\n";
    return 1;
  }

  for(int i = 0; i < slacks.size(); i++){
    if(skip_start <= node_count + i && node_count + i <= skip_end)
      continue;
    if(fabs(slacks[i]) >= LP::EPSILON){
      delset[node_count + i] = 1;
      num_removed++;
    }
  }

  rval = PSEPlp_delsetrows(&m_lp, &delset[0]);
  if(rval){
    cerr << "Problem in LPPrune::prune_with_skip with delsetrows\n";
    return 1;
  }

  return 0;
}
