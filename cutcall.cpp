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
  dominos.test_build_collection();
  
  return 0;
}

int PSEP_SimpleDP::in_subtour_poly(bool *result_p){
  int ecount = support_ecap.size(), ncount = support_ecap.size() / 2;
  int end0 = 0;
  double cutval = 2;
  *result_p = false;
  
  for(int end1 = 1; end1 < ncount; end1++){
    if(CCcut_mincut_st(ncount, ecount, &support_elist[0], &support_ecap[0],
		       end0, end1, &cutval, (int **) NULL, (int *) NULL)){
      cerr << "Problem in SimpleDP::separate with Concorde st-cut" << endl;
      return 1;
    }

    if(cutval < 2)
      return 0;
  }

  *result_p = true;
  return 0;
}
