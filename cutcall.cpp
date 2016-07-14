#include "cutcall.h"

using namespace std;

int PSEP_Cutcall::segment(int *num_added_p){
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
    G_Utils::get_delta(segnodes, edges, &deltacount, delta,
		       edge_marks);

    rval = segments.add_cut(deltacount, delta);
    if(rval) goto CLEANUP;
      
    (*num_added_p)++;
  }

  if(*num_added_p == 0) rval = 2;
    
 CLEANUP:
  if(rval == 1)
    cerr << "Entry point: PSEP_Cutcall::segment" << endl;
  return rval;
}

int PSEP_Cutcall::blossom(const int max_cutcount, int *num_added_p){
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
    G_Utils::get_delta(hnodes, edges, &deltacount, delta, edge_marks);

    rval = blossoms.add_cut(deltacount, delta, cutedge);
    if(rval) goto CLEANUP;
    (*num_added_p)++;
  }

  if(*num_added_p == 0)
    rval = 2;


 CLEANUP:
  if(rval == 1)
    cerr << "Error entry point: blossom_cutcall" << endl;
  return rval;
}

int PSEP_Cutcall::simpleDP(const int max_cutcount, int *num_added_p){
  bool in_sep = false;

  if(dominos.in_subtour_poly(&in_sep))
    return 1;
  
  if(in_sep){
    cout << "Solution is in the subtour polytope, building collection..\n";
    dominos.test_build_collection();
  } else
    cout << "Solution is not in subtour polytope, exiting\n";

  return 0;
}
