#include "segcuts.h"

using namespace std;

void PSEP_Segcuts::separate(){
  int ncount = G_s.node_count;
  int current_start, current_end, current_size;
  SNode current_snode;
  double lhs;

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
      
      if(lhs > current_size - 1 && (fabs(lhs - (current_size - 1)) >= 0.002)){
	seg new_seg(i, j, lhs - (current_size - 1));
	if(pq.size() == 250){
	  if(new_seg < pq.top()){
	    pq.pop();
	    pq.push(new_seg);
	  } 
	} else
	  pq.push(new_seg);
      }
    }
    for(int l = i; l <= j; l++)
      edge_marks[best_tour_nodes[l]] = 0;
  }
}

int PSEP_Segcuts::add_cut(const int deltacount, vector<int> &delta){
  int rval = 0, newrows = 1, newnz = deltacount;
  int rmatbeg[1] = {0};
  char sense[1] = {'G'};
  double rhs[1] = {2.0};
  vector<double> rmatval(deltacount, 1.0);

  rval = PSEPlp_addrows (&m_lp, newrows, newnz, rhs, sense, rmatbeg,
			 &delta[0], &rmatval[0]);

  if(rval)
    cerr << "Entry point: PSEP_Segments::add_cut" << endl;
  return rval;
}
