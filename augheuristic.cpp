#include "augheuristic.h"

int PSEP_AugHeuristic::add_clamp(){
  int rval = 0;
  double min_dif = 1.0;
  int best_edge = -1;
  double newbound;

  for(int i = 0; i < support_indices.size(); i++){
    if(support_ecap[i] > 1 - LP::EPSILON) continue;
    int index = support_indices[i];
    if(fabs(support_ecap[i] - best_tour_edges[index]) < min_dif){
      min_dif = fabs(support_ecap[i] - best_tour_edges[index]);
      best_edge = index;
    }
  }

  switch(best_tour_edges[best_edge]){
  case 0:
    newbound = 0.0;
    break;
  case 1:
    newbound = 1.0;
  }

  rval = PSEPlp_setbnd(&m_lp, best_edge, 'B', newbound);
  if(rval)
    std::cerr << "Error entry point: AugHeuristic::add_clamps()\n";
  else{
    std::cout << "Added heuristic clamp on edge "<< best_edge
	      << ", clamped to " << newbound << " ("
	      << best_tour_edges[best_edge] << "), dif: "
	      << min_dif << "\n";
    clamp_edges.push_back(best_edge);
  }

  return rval;
}

int PSEP_AugHeuristic::clear_clamps(){
  int rval = 0;
  while(active()){
    int current_edge = clamp_edges.back();
    clamp_edges.pop_back();
    
    rval = PSEPlp_setbnd(&m_lp, current_edge, 'L', 0.0);
    if(rval) goto CLEANUP;
    
    rval = PSEPlp_setbnd(&m_lp, current_edge, 'U', 1.0);
    if(rval) goto CLEANUP;
  }

 CLEANUP:
  if(rval){
    std::cerr << "Error entry point: AugHeuristic::clear_clamps()\n";
    clamp_edges.clear();
  } else
    std::cout << "Removed all clamps successfully.\n";
  return rval;
}
