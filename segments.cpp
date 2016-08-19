#include<iostream>

#include "segments.h"

using namespace std;
using namespace PSEP;

int Cut<seg>::separate(){
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
	if(!best || (lhs - (current_size - 1)) > best->viol)
	  best.reset(new seg(i, j,lhs - (current_size - 1)));
      }
    }
    for(int l = i; l <= j; l++)
      edge_marks[best_tour_nodes[l]] = 0;
  }

  if(!best)
    return 2;
  
  return 0;
}

int Cut<seg>::parse_coeffs(){
  if(!best){
    cerr << "Cuts<seg>::parse_coeffs tried to parse null pointer\n";
    return 1;
  }

  deltacount = 0;

  std::vector<int> segnodes;

  for(int i = best->start; i <= best->end; i++)
    segnodes.push_back(best_tour_nodes[i]);

  GraphUtils::get_delta(segnodes, edges, &deltacount, delta, edge_marks);

  if(deltacount == 0){
    cerr << "Cuts<seg>::parse_coeffs returned no edges\n";
    return 1;
  }
  
  return 0;
}

int Cut<seg>::add_cut(){
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

int Cut<seg>::cutcall(){
  int rval = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = parse_coeffs();
  if(rval) goto CLEANUP;

  rval = add_cut();
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cuts<seg>::cutcall()\n";
  best.reset(NULL);
  return rval;
}
