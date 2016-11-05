#include "blossoms.hpp"

#include <iostream>

using std::vector;
using std::cerr;

#define PSEP_OMP_2MATCH

#ifdef PSEP_OMP_2MATCH
#warning USING OMP VERSION
#else
#warning USING SERIAL VERSION
#endif

namespace PSEP {

#ifdef PSEP_OMP_2MATCH

int Cut<blossom>::separate(){
  int rval = 0;
  int ncount = m_graph.node_count;

  try { cut_ecap.resize(support_ecap.size()); } catch(...){
    rval = 1; PSEP_GOTO_CLEANUP("Couldn't resize cut_ecap, ");
  }

  for(int i = 0; i < support_indices.size(); i++){
    if(best_tour_edges[support_indices[i]] == 0)
      cut_ecap[i] = support_ecap[i];
    else
      cut_ecap[i] = 1 - support_ecap[i];
  }

  #pragma omp parallel for
  for(int i = 0; i < support_indices.size(); i++){
    int cut_edge_index = support_indices[i];
    int best_tour_entry = best_tour_edges[cut_edge_index];
    int end0 = support_elist[2 * i];
    int end1 = support_elist[(2 * i) + 1];
    
    int cutcount = 0;
    int *cut_nodes = (int *) NULL;
    double orig_weight = 1.0, changed_weight = 1.0, cutval = 1.0;

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

    rval = CCcut_mincut_st(ncount, support_indices.size(),
			   &support_elist[0], &cut_ecap[0], end0, end1,
			   &cutval, &cut_nodes, &cutcount);
    if(rval){
      cerr << "Problem in blossom::separate w st-cut\n";
      i = support_indices.size();
      cutval = 1.0;
    }

    if(cutval < 1 - LP::EPSILON && cutcount >= 3 &&
       cutcount <= (ncount - 3)){
      
      vector<int> handle;
      for(int j = 0; j < cutcount; j++){
	handle.push_back(cut_nodes[j]);
      }

      blossom new_cut(handle, cut_edge_index, cutval);

      #pragma omp critical
      {
      try { //if it is a better cut it goes to the front for immediate adding
	if(local_q.empty() || cutval <= local_q.peek_front().cut_val){
	  local_q.push_front(new_cut);
	}
	else { //it goes to the back for use in the pool if applicable
	  local_q.push_back(new_cut);
	}
      } catch (...) {
	rval = 1; //PSEP_GOTO_CLEANUP("Problem pushing new cut to queue, ");
      }
      }
      
    }

    if(rval){
      cerr << "Problem pushing new cut to queue. ";
      i = support_indices.size();
    }

    cut_ecap[i] = orig_weight;
    CC_IFFREE(cut_nodes, int);
  }

  if(local_q.empty()) rval = 2;


 CLEANUP:
  return rval;
}

#else

int Cut<blossom>::separate(){
  int rval = 0;
  int ncount = m_graph.node_count;

  try { cut_ecap.resize(support_ecap.size()); } catch(...){
    rval = 1; PSEP_GOTO_CLEANUP("Couldn't resize cut_ecap, ");
  }

  for(int i = 0; i < support_indices.size(); i++){
    if(best_tour_edges[support_indices[i]] == 0)
      cut_ecap[i] = support_ecap[i];
    else
      cut_ecap[i] = 1 - support_ecap[i];
  }

  for(int i = 0; i < support_indices.size(); i++){
    int cut_edge_index = support_indices[i];
    int best_tour_entry = best_tour_edges[cut_edge_index];
    int end0 = support_elist[2 * i];
    int end1 = support_elist[(2 * i) + 1];
    
    int cutcount = 0;
    int *cut_nodes = (int *) NULL;
    double orig_weight = 1.0, changed_weight = 1.0, cutval = 1.0;

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

    rval = CCcut_mincut_st(ncount, support_indices.size(),
			   &support_elist[0], &cut_ecap[0], end0, end1,
			   &cutval, &cut_nodes, &cutcount);
    if(rval){
      cerr << "Problem in blossom::separate w st-cut\n";
      goto CLEANUP;
    }

    if(cutval < 1 - LP::EPSILON && cutcount >= 3 &&
       cutcount <= (ncount - 3)){
      
      vector<int> handle;
      for(int j = 0; j < cutcount; j++){
	handle.push_back(cut_nodes[j]);
      }

      blossom new_cut(handle, cut_edge_index, cutval);

      try { //if it is a better cut it goes to the front for immediate adding
	if(local_q.empty() || cutval <= local_q.peek_front().cut_val){
	  local_q.push_front(new_cut);
	}
	else { //it goes to the back for use in the pool if applicable
	  local_q.push_back(new_cut);
	}
      } catch (...) {
	rval = 1; PSEP_GOTO_CLEANUP("Problem pushing new cut to queue, ");
      }      
    }

    cut_ecap[i] = orig_weight;
    CC_IFFREE(cut_nodes, int);
  }

  if(local_q.empty()) rval = 2;


 CLEANUP:
  return rval;
}

#endif

int Cut<blossom>::build_hypergraph(const blossom &blossom_cut){
  int rval = 0;
  int cutedge = blossom_cut.cut_edge;
  int num_teeth = 0;
  deltacount = 0;
  vector<vector<int>> node_sets;
  vector<Edge> &edges(m_graph.edges);

  
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

  num_teeth = node_sets.size() - 1;

  if(num_teeth %2 == 1){
    try {
      HyperGraph newcut(node_sets, HyperGraph::CutType::Blossom);
      external_q.push_back(newcut);
    } catch (...) {
      rval = 1; PSEP_GOTO_CLEANUP("Couldn't push to external queue, ");
    }
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

  if(external_q.empty()) rval = 2;

 CLEANUP:
  if(rval == 1)
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
  return rval;
}

}
