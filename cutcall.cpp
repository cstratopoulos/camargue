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
    cerr << "Error entry point: PSEP_cutcall::blossom" << endl;
  return rval;
}

int PSEP_Cutcall::simpleDP(const int max_cutcount, int *num_added_p){
  int rval = 0;
  *num_added_p = 0;
  vector<double> agg_coeffs(m_lp_edges.size(), 0);
  // double rhs, lhs;
  // int RHS;
  int ncount, ecount;
  vector<int> cut_node_marks;
  vector<int> domino_delta;
  int deltacount = 0;

  if(max_cutcount == 0){
    rval = 2;
    goto CLEANUP;
  }

  cout << "simpleDP is calling dominos.separate!\n";
  rval = dominos.separate(max_cutcount);
  if(rval) goto CLEANUP;

  if(dominos.toothlists.empty()){
    rval = 2; goto CLEANUP;
  }

  cout << "dominos.toothlists.size is " << dominos.toothlists.size() << "\n";

  ncount = dominos.light_nodes.size();
  ecount = dominos.cut_ecap.size();
  domino_delta.resize(ecount);
  cut_node_marks.resize(ncount);

  for(int i = 0; i < dominos.toothlists.size(); i++){
    double rhs = 0.0;
    double lhs = 0.0;
    int RHS = 0;
    vector<int> &current_toothlist = dominos.toothlists[i];
    for(int j = 0; j < agg_coeffs.size(); j++)
      agg_coeffs[j] = 0;
    
    G_Utils::get_delta(current_toothlist.size(),
		       &current_toothlist[0], ecount,
		       &(dominos.cut_elist)[0], &deltacount, &domino_delta[0],
		       &cut_node_marks[0]);

    cout << "============================\n";
    cout << "NOW PARSING TOOTHLIST " << i << "\n";

    dominos.parse_domino(deltacount, domino_delta, agg_coeffs, &rhs);

    for(int j = 0; j < agg_coeffs.size(); j++)
      agg_coeffs[j] /= 2;

    rhs /= 2;
    RHS = (int) rhs;
    rhs = RHS;

    for(int j = 0; j < m_lp_edges.size(); j++)
      lhs += m_lp_edges[j] * agg_coeffs[j];

    if(lhs <= rhs){
      cerr << "Bad inequality won't be added, lhs: " << lhs << ", rhs: "
	   << rhs << "\n";
      continue;
    } else {
      cout << "Found simple DP with lhs: "
	   << lhs << ", rhs: " << rhs << "......";
    }
    

    rval = dominos.add_cut(agg_coeffs, rhs);
    if(rval)
      goto CLEANUP;
    cout << "Added inequality.\n";

    (*num_added_p)++;
  }


 CLEANUP:
  if(rval == 1)
    cerr << "Error entry point: PSEP_cutcall::simpleDP" << endl;
  return 0;
}

int PSEP_Cutcall::in_subtour_poly(bool *result_p){
  int ecount = support_ecap.size(), ncount = best_tour_nodes.size();  
  int end0 = 0;
  double cutval = 2;
  *result_p = false;
  
  for(int end1 = 1; end1 < ncount; end1++){
    if(CCcut_mincut_st(ncount, ecount, &support_elist[0], &support_ecap[0],
		       end0, end1, &cutval, (int **) NULL, (int *) NULL)){
      cerr << "Problem in Cutcall::in_subtour_poly w Concorde st-cut" << endl;
      return 1;
    }

    if(cutval < 2)
      return 0;
  }

  *result_p = true;
  return 0;
}
