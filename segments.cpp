#include<iostream>

#include "segments.hpp"

using namespace std;
using namespace PSEP;

int Cut<seg>::separate(){
  int ncount = G_s.node_count;
  int current_start, current_end, current_size;
  SNode current_snode;
  double lhs, best_viol = 0;

  for(int i = 0; i < ncount - 2; i++){
    current_start = best_tour_nodes[i];
    edge_marks[current_start] = 1;
    current_size = 1;
    lhs = 0;
    int j;
    
    for(j = i + 1; (j < ncount - 1 && (++current_size) <= ncount / 2); j++){
      current_end = best_tour_nodes[j];
      edge_marks[current_end] = 1;
      current_snode = G_s.nodelist[current_end];
      
      for(int k = 0; k < current_snode.s_degree; k++)
	if(edge_marks[current_snode.adj_objs[k].other_end] == 1)
	  lhs += current_snode.adj_objs[k].lp_weight;

      double viol = fabs(lhs - (current_size - 1));
      if(lhs > current_size - 1 && viol >= 0.002){
	if(viol >= best_viol){
	  best_viol = viol;
	  seg newbest(i, j, viol);
	  try { local_q.push_front(newbest); } catch(...) {
	    cerr << "Problem pushing cut to local queue in Cut<seg>separate\n";
	    return 1;
	  }
	}
      }
    }
    for(int l = i; l <= j; l++)
      edge_marks[best_tour_nodes[l]] = 0;
  }

  if(local_q.empty())
    return 2;
  
  return 0;
}

int Cut<seg>::build_hypergraph(const seg& seg_cut){
  int rval = 0;
  
  if(HyperGraph::same_tour(best_tour_nodes)){
    try {
    HyperGraph newcut(IntPair(seg_cut.start, seg_cut.end));
    external_q.push_front(newcut);
    } catch (...) {
      rval = 1; PSEP_GOTO_CLEANUP("Problem adding hypergraph to queue, ");
    }
  } else {
    vector<vector<int>> segment_nodes;
    segment_nodes.resize(1);
    try {
      for(int i = seg_cut.start; i <= seg_cut.end; i++)
	segment_nodes[0].push_back(best_tour_nodes[i]);

      HyperGraph newcut(segment_nodes, HyperGraph::CutType::Segment);
      external_q.push_front(newcut);
    } catch (...){
      rval = 1; PSEP_GOTO_CLEANUP("Problem adding hypergraph to queue, ");
    }
  }

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<seg>::build_hypergraph\n";
  return rval;
}

int Cut<seg>::add_cuts(){
  int rval = 0;
  
  while(!local_q.empty()){
    rval = build_hypergraph(local_q.peek_front());
    if(rval) goto CLEANUP;
    local_q.pop_front();
  }

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<seg>::add_cuts\n";
  return rval;
}

int Cut<seg>::cutcall(){
  int rval = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = add_cuts();
  if(rval) goto CLEANUP;
  
 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cuts<seg>::cutcall()\n";
  return rval;
}
