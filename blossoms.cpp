#include<iostream>

#include "blossoms.hpp"

using namespace std;
using namespace PSEP;

int Cut<blossom>::separate(){
  int rval = 0;
  int cut_edge_index;
  int end0, end1;
  int best_tour_entry;
  int ncount = 2 * support_indices.size();
  double orig_weight, changed_weight, cutval, min_cutval = 1.0;
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
      if(cutval <= min_cutval){
	min_cutval = cutval;
	vector<int> handle;
	for(int j = 0; j < cutcount; j++){
	  handle.push_back(cut_nodes[j]);
	}

	blossom new_cut(handle, cut_edge_index, cutval);
	try { local_q.push_front(new_cut); } catch (...) {
	  rval = 1; PSEP_GOTO_CLEANUP("Problem pushing new cut to queue, ");
	}
      }
    }

    cut_ecap[i] = orig_weight;
    CC_IFFREE(cut_nodes, int);
  }

  if(local_q.empty()) rval = 2;

 CLEANUP:
  CC_IFFREE(cut_nodes, int);
  return rval;
}


int Cut<blossom>::build_hypergraph(const blossom &blossom_cut){
  int rval = 0;
  int cutedge = blossom_cut.cut_edge;
  deltacount = 0;
  vector<vector<int>> node_sets;
  
  try { node_sets.push_back(blossom_cut.handle); } catch (...) {
    rval = 1; PSEP_GOTO_CLEANUP("Problem pushing back node set, ");
  }

  GraphUtils::get_delta(blossom_cut.handle, edges, &deltacount, delta,
			edge_marks);

  if(deltacount == 0){
    rval = 1; PSEP_GOTO_CLEANUP("Get delta returned zero handle edges, ");
  }

  switch(best_tour_edges[cutedge]){
  case 0:
    for(int i = 0; i < deltacount; i++)
      if(m_lp_edges[delta[i]] > LP::EPSILON &&
	 (best_tour_edges[delta[i]] == 1 || delta[i] == cutedge)){
	int edge_index = delta[i];
	vector<int> new_tooth{edges[edge_index].end[0],
	    edges[edge_index].end[1]};
	try { node_sets.push_back(new_tooth); } catch (...) {
	  rval = 1; PSEP_GOTO_CLEANUP("Problem pushing back node set, ");
	}
      }
    break;

  case 1:
    for(int i = 0; i < deltacount; i++)
      if(m_lp_edges[delta[i]] > LP::EPSILON &&
	 (best_tour_edges[delta[i]] == 1 && delta[i] != cutedge)){
	int edge_index = delta[i];
	vector<int> new_tooth{edges[edge_index].end[0],
	    edges[edge_index].end[1]};
	try { node_sets.push_back(new_tooth); } catch (...) {
	  rval = 1; PSEP_GOTO_CLEANUP("Problem pushing back node set, ");
	}
      }
  }

  try {
  HyperGraph newcut(node_sets, HyperGraph::CutType::Blossom);
  external_q.push_front(newcut);
  } catch (...) {
    rval = 1; PSEP_GOTO_CLEANUP("Couldn't push to external queue, ");
  }

  

 CLEANUP:
  if(rval)
    cerr << "Cut<blossom>::build_hypergraph failed\n";
  return rval;
}

int Cut<blossom>::add_cuts(){
  int rval = 0;

  while(!local_q.empty()){
    rval = build_hypergraph(local_q.peek_front());
    if(rval) goto CLEANUP;
    local_q.pop_front();
  }


 CLEANUP:
  if(rval)
    cerr << "problem in Cut<blossom>::add_cuts\n";
  return rval;
}

int Cut<blossom>::cutcall(){
  int rval = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = add_cuts();
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cuts<blossom>::cutcall\n";
  best.reset(NULL);
  return rval;
}
