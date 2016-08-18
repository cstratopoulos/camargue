#include<iostream>

#include "blossoms.h"

using namespace std;
using namespace PSEP;

int Cut<blossom>::separate(){
  int rval = 0;
  int cut_edge_index;
  int end0, end1;
  int best_tour_entry;
  int ncount = 2 * support_indices.size();
  double orig_weight, changed_weight, cutval;
  int *cut_nodes = (int *) NULL;
  int cutcount = 0;

  cut_ecap.resize(support_ecap.size());

  for(int i = 0; i < support_indices.size(); i++){
    if(best_tour_edges[support_indices[i]] == 0)
      cut_ecap[i] = support_ecap[i];
    else
      cut_ecap[i] = 1 - support_ecap[i];
  }

  for(int i = 0; i < support_indices.size(); i++){
    cut_edge_index = support_indices[i];
    best_tour_entry = best_tour_edges[cut_edge_index];

    switch(best_tour_entry){
    case(0):
      orig_weight = support_ecap[i];
      changed_weight = 1 - support_ecap[i];
      break;
    case(1):
      orig_weight = 1 - support_ecap[i];
      changed_weight = support_ecap[i];
    }

    cut_ecap[i] = changed_weight;

    end0 = support_elist[2 * i];
    end1 = support_elist[(2 * i) + 1];

    rval = CCcut_mincut_st(ncount, support_indices.size(),
			   &support_elist[0], &cut_ecap[0], end0, end1,
			   &cutval, &cut_nodes, &cutcount);
    if(rval){
      cerr << "Problem in blossom::separate w st-cut" << endl;
      goto CLEANUP;
    }

    if(cutval < 1 - LP::EPSILON){
      if(!best || cutval < best->cut_val){
	vector<int> handle;
	for(int j = 0; j < cutcount; j++){
	  handle.push_back(cut_nodes[j]);
	}

	best.reset(new blossom(handle, cut_edge_index, cutval));
      }
    }

    cut_ecap[i] = orig_weight;
  }

  if(!best) rval = 2;

 CLEANUP:
  CC_IFFREE(cut_nodes, int);
  return rval;
}

int Cut<blossom>::parse_coeffs(){
  if(!best){
    cerr << "Cuts<blossom>::parse_coeffs tried to parse null pointer\n";
    return 1;
  }

  deltacount = 0;

  G_Utils::get_delta(best->handle, edges, &deltacount, delta, edge_marks);

  if(deltacount == 0){
    cerr << "Cuts<blossom>::parse_coeffs returned zero edges\n";
    return 1;
  }

  return 0;
}

int Cut<blossom>::add_cut(){
  int rval = 0, newrows = 1, newnz = deltacount;
  int rmatbeg[1] = {0};
  char sense[1] = {'G'};
  double rhs[1];
  int num_teeth = 0;
  vector<double> rmatval(deltacount, 1.0);
  int cutedge = best->cut_edge;

  switch(best_tour_edges[cutedge]){
  case 0:
    for(int i = 0; i < deltacount; i++)
      if(m_lp_edges[delta[i]] > LP::EPSILON &&
	 (best_tour_edges[delta[i]] == 1 || delta[i] == cutedge)){
	rmatval[i] = -1.0;
	num_teeth++;
      }
    break;
  case 1:
    for(int i = 0; i < deltacount; i++)
      if(m_lp_edges[delta[i]] > LP::EPSILON &&
	 (best_tour_edges[delta[i]] == 1 && delta[i] != cutedge)){
	rmatval[i] = -1.0;
	num_teeth++;
      }
  }

  rhs[0] = 1 - num_teeth;  

  rval = PSEPlp_addrows (&m_lp, newrows, newnz, rhs, sense, rmatbeg,
			 &delta[0], &rmatval[0]);

  if(rval)
    cerr << "Entry point: Cut<blossom>::add_cut" << endl;
  return rval;
}

int Cut<blossom>::cutcall(){
  int rval = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = parse_coeffs();
  if(rval) goto CLEANUP;

  rval = add_cut();
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cuts<blossom>::cutcall\n";
  best.reset(NULL);
  return rval;
}
