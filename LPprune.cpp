#include<iostream>
#include<algorithm>

#include "LPprune.hpp"

using namespace std;
using namespace PSEP::LP;

constexpr int MAX_AGE = 100;
constexpr double EPSILON_PRUNE = 0.000001;

int CutPrune::aggressive_prune(){
  int numrows = PSEPlp_numrows(&m_lp);
  
  dual_vars.resize(PSEPlp_numrows(&m_lp) - node_count, 1);
  ages.resize(PSEPlp_numrows(&m_lp) - node_count, 0);
  delset.resize(numrows, 0);

  if(ages.empty()) return 0;
  
  bool found_old = false;  
  int i;

  if(dual_measure()){
    cerr << "Problem in CutPrune::aggressive_prune\n";
    return 1;
  }
  
  for(i = 0; i < ages.size(); i++)
    if(ages[i] >= MAX_AGE){
      found_old = true;
      break;
    }
  if(!found_old) return 0;

  
  cout << "Gonna remove row:";
  for( ; i < ages.size(); i++)
    if(ages[i] >= MAX_AGE){
      delset[i + node_count] = 1;
      cout << (i + node_count) << ", age " << ages[i] << "\n";
    }

  if(PSEPlp_delsetrows(&m_lp, &delset[0])){
    cerr << "Problem in CutPrune::aggressive_prune\n";
    return 1;
  }

  ages.erase(remove_if(ages.begin(), ages.end(),
		       [](int n){ return n >= MAX_AGE; }),
	     ages.end());
  dual_vars.resize(ages.size());

  for(int j = node_count; j < delset.size(); j++)
    delset[j] = 0;

  cout << "Aggressive removed "
       << (numrows - PSEPlp_numrows(&m_lp)) << " cuts from LP\n";
  return 0;
}

int CutPrune::prune_cuts(int &num_removed){
  int rval = 0;
  num_removed = 0;
  int rowcount = PSEPlp_numrows(&m_lp);
  vector<int> delset(rowcount, 0);
  vector<double> slacks(rowcount - node_count, 0);

  rval = PSEPlp_getslack(&m_lp, &slacks[0], node_count, rowcount - 1);
  if(rval){
    cerr << "Problem in CutPrune::prune_cuts with getslack\n";
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
    cerr << "Problem in CutPrune::prune_cuts with delsetrows\n";
    return 1;
  }

  return 0;
}

int CutPrune::prune_with_skip(int &num_removed, IntPair skiprange,
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
    cerr << "Problem in CutPrune::prune_with_skip with getslack\n";
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
    cerr << "Problem in CutPrune::prune_with_skip with delsetrows\n";
    return 1;
  }

  return 0;
}

int CutPrune::dual_measure(){
  if(PSEPlp_getpi(&m_lp, &dual_vars[0], node_count,
		  PSEPlp_numrows(&m_lp) -1)){
    cerr << "Problem in CutPrune::dual_measure\n";
    return 1;
  }

  for(int i = 0; i < dual_vars.size(); i++){
    if(dual_vars[i] < EPSILON_PRUNE)
      ages[i]++;
    else
      ages[i] = 0;
  }

  return 0;
}
