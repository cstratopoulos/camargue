#include "segments.hpp"
#include "PSEP_util.hpp"

#include<iostream>

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

using std::vector;
using std::cout;
using std::cerr;

namespace PSEP {

int Cut<seg>::linsub_callback(double cut_val, int cut_start, int cut_end,
			      void *cb_data)
{
  Cut<seg> *this_p = (Cut<seg> *) cb_data;
  CutQueue<seg> &lq = this_p->local_q;

  if(lq.empty() || cut_val <= lq.peek_front().cutval){
    try{ lq.push_front(seg(cut_start, cut_end, cut_val)); } catch(...){
      cerr << "Couldn't push back new seg\n";
      return 1;
    }
  }

  return 0;
}

int Cut<seg>::separate()
{
  int rval = 0;
  vector<int> endmark;

  vector<int> cut_elist;

  try { cut_elist.resize(support_elist.size()); } catch(std::bad_alloc &) {
    PSEP_SET_GOTO(rval, "Out of space for cut elist. ");
  }

  for(int i = 0; i < cut_elist.size(); i++)
    cut_elist[i] = perm[support_elist[i]];


  try { endmark = vector<int>(G_s.node_count, CC_LINSUB_BOTH_END); }
  catch(...){
    rval = 1; PSEP_GOTO_CLEANUP("Couldn't allocate endmark. ");
  }

  rval = CCcut_linsub(G_s.node_count, G_s.edge_count, &endmark[0],
		      &cut_elist[0], &support_ecap[0], 2.0 - Epsilon::Cut,
		      (void *) this, linsub_callback);
  PSEP_CHECK_RVAL(rval, "CCcut_linsub failed. ");

  if(local_q.empty())
    rval = 2;
  

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<seg>::separate\n";
  return rval;
}

int Cut<seg>::build_hypergraph(const seg& seg_cut){
  int rval = 0;
  
  if(HyperGraph::same_tour(best_tour_nodes)){
    try {
    HyperGraph newcut(IntPair(seg_cut.start, seg_cut.end));
    external_q.push_back(newcut);
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
      external_q.push_back(newcut);
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

}
