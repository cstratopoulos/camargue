#include<vector>
#include<iostream>

#include "segments2.h"
#include "Graph.h"
#include "lp.h"

using namespace std;
using namespace PSEP;

template<>
int Cuts<seg>::separate(){
  SupportGraph &G_s = SupportGroup.G_s;
  std::vector<int> &edge_marks = GraphGroup.edge_marks;
  std::vector<int> &best_tour_nodes = BestGroup.best_tour_nodes;

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
	if((lhs - (current_size - 1)) > best->viol)
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

template<>
int Cuts<seg>::get_coefficients(int *deltacount){
  if(!best){
    cerr << "Cuts<seg>::get_coefficients tried to parse null pointer\n";
    return 1;
  }

  *deltacount = 0;

  std::vector<int> &best_tour_nodes = BestGroup.best_tour_nodes;
  std::vector<Edge> &edges = GraphGroup.m_graph.edges;
  std::vector<int> &delta = GraphGroup.delta;
  std::vector<int> &edge_marks = GraphGroup.edge_marks;

  std::vector<int> segnodes;

  for(int i = best->start; i <= best->end; i++)
    segnodes.push_back(best_tour_nodes[i]);

  G_Utils::get_delta(segnodes, edges, deltacount, delta, edge_marks);
  
  return 0;
}

template<>
int Cuts<seg>::add_cut(const int deltacount){
  PSEPlp &m_lp = LPGroup.m_lp;
  
  std::vector<int> &delta = GraphGroup.delta;
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

template<>
int Cuts<seg>::cut_call(){
  int rval = 0, deltacount = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = get_coefficients(&deltacount);
  if(rval || deltacount < 0) goto CLEANUP;

  rval = add_cut(deltacount);
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval)
    cerr << "Problem in Cuts<seg>::cutcall()\n";
  best.reset(NULL);
  return rval;
}
