#include "tsp_solver.h"

#include<iomanip>

using namespace std;

int TSP_Solver::seg_cutcall(int *num_added_p){
  int rval = 0;
  *num_added_p = 0;

  segments.separate();

  while(!segments.q_empty()){
    int start, end; double viol;
    segments.pop(&start, &end, &viol);

    vector<int> segnodes;
    for(int i = start; i <= end; i++)
      segnodes.push_back(best_tour_nodes[i]);
    
    int deltacount;
    G_Utils::get_delta(segnodes, m_graph.edges, &deltacount, delta,
		       edge_marks);

    rval = segments.add_cut(deltacount, delta);
    if(rval) goto CLEANUP;
      
    (*num_added_p)++;
  }

  if(*num_added_p == 0) rval = 2;
  else
    cout << "Added " << *num_added_p << " segment cuts" << endl;
    
 CLEANUP:
  if(rval == 1)
    cerr << "Entry point: TSP_Solver::seg_cutcall" << endl;
  return rval;
}

int TSP_Solver::blossom_cutcall(const int max_cutcount, int *num_added_p){
  int rval = 0;
  *num_added_p = 0;

  if(max_cutcount == 0)
    goto CLEANUP;
  
  rval = blossoms.separate(max_cutcount);
  if(rval) goto CLEANUP;

  while(!blossoms.q_empty()){
    vector<int> hnodes;
    int cutedge;
    double cutval;
    blossoms.pop(hnodes, &cutedge, &cutval);

    int deltacount;
    G_Utils::get_delta(hnodes, m_graph.edges, &deltacount, delta, edge_marks);

    rval = blossoms.add_cut(deltacount, delta, cutedge);
    if(rval) goto CLEANUP;
    (*num_added_p)++;
  }

  if(*num_added_p == 0){
    rval = 2;
  } else
    cout << "Added " << *num_added_p << " blossoms in total" << endl;

 CLEANUP:
  if(rval == 1)
    cerr << "Error entry point: blossom_cutcall" << endl;
  return rval;
}
