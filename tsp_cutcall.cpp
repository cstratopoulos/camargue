#include "tsp_solver.h"

#include<iomanip>

using namespace std;

int TSP_Solver::seg_cutcall(int *num_added_p){
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
    G_Utils::get_delta(segnodes, m_graph.edges, &deltacount, delta,
		       edge_marks);

    rval = segments.add_cut(deltacount, delta);
    if(rval) goto CLEANUP;
      
    (*num_added_p)++;
  }

  if(*num_added_p == 0) rval = 2;
  else
    cout << "Added " << *num_added_p << " segment cuts" << endl;
    
 CLEANUP:
  if(rval == 1)
    cerr << "Entry point: TSP_Solver::seg_cutcall" << endl;
  return rval;
}

int TSP_Solver::blossom_cutcall(const int max_cutcount, int *num_added_p){
  int rval = 0;
  *num_added_p = 0;
  
  rval = blossoms.separate(max_cutcount);
  if(rval) goto CLEANUP;

  while(!blossoms.q_empty()){
    cout << ">>>>>> POPPING A BLOSSOM <<<<<<<" << endl;
    vector<int> hnodes;
    int cutedge;
    double cutval;
    blossoms.pop(hnodes, &cutedge, &cutval);

    cout << "Considering blossom with cutval "
	 << setprecision(10) << cutval << endl;
    cout << "Cut edge is " << cutedge;
    cout << ", is cutedge in best tour: " << best_tour_edges[cutedge] << endl;

    int deltacount;
    G_Utils::get_delta(hnodes, m_graph.edges, &deltacount, delta, edge_marks);

    cout << deltacount << " edges in delta, nonzero ones: " << endl;
    for(int i = 0; i < deltacount; i++){
      if(m_lp_edges[delta[i]] > 0 || best_tour_edges[delta[i]] > 0){
	cout << "LP val: " << m_lp_edges[delta[i]] << ", ";
	cout << "In best tour: " << best_tour_edges[delta[i]] << ", ";
	cout << "Cut ecap: ";
	if(delta[i] == cutedge){
	  if(best_tour_edges[delta[i]] == 1)
	    cout << support_ecap[i] << ", ";
	  else
	    cout << 1 - support_ecap[i] << ", ";
	} else {
	  if(best_tour_edges[i] == 1)
	    cout <<  1 - support_ecap[i] << ", ";
	  else
	    cout << support_ecap[i] << ", ";
	}
	m_graph.print_edge(delta[i]);
      }

    }
    rval = blossoms.add_cut(deltacount, delta, cutedge);
    if(rval) goto CLEANUP;
    (*num_added_p)++;
  }

  if(*num_added_p == 0){
    rval = 2;
  } else
    cout << "Added " << *num_added_p << " blossoms in total" << endl;

 CLEANUP:
  if(rval == 1)
    cerr << "Error entry point: blossom_cutcall" << endl;
  return rval;
}
