#include "augheuristic.h"

using namespace std;

int PSEP_AugHeuristic::add_clamp(){
  int rval = 0;
  double min_dif = 1.0;
  int best_edge = -1, found_ind = -1;
  double newbound;

  for(int i = 0; i < support_indices.size(); i++){
    if(support_ecap[i] > 1 - LP::EPSILON) continue;
    int index = support_indices[i];
    if(fabs(support_ecap[i] - best_tour_edges[index]) < min_dif){
      min_dif = fabs(support_ecap[i] - best_tour_edges[index]);
      best_edge = index;
      found_ind = i;
    }
  }

  switch(best_tour_edges[best_edge]){
  case 0:
    newbound = 0.0;
    break;
  case 1:
    newbound = 1.0;
  }

  rval = PSEPlp_clampbnd(&m_lp, best_edge, newbound);
  if(rval)
    std::cerr << "Error entry point: AugHeuristic::add_clamps()\n";
  else{
    std::cout << "||||| Added heuristic clamp on edge "<< best_edge
	      << ", clamped to " << newbound << " ("
	      << best_tour_edges[best_edge] << "), dif: "
	      << min_dif << "\n";
    std::cout << "Check support index: " << support_indices[found_ind]
	      << ", cap: " << support_ecap[found_ind] << "\n";
    clamp_edges.push_back(best_edge);
  }

  return rval;
}

int PSEP_AugHeuristic::clear_clamps(){
  int rval = 0;
  std::vector<double> bounds(clamp_edges.size(), 0.0);
  std::vector<char> lower_or_upper(clamp_edges.size(), 'L');

  rval = PSEPlp_relaxbds(&m_lp, clamp_edges.size(), &clamp_edges[0],
			 &lower_or_upper[0], &bounds[0]);
  if(rval) goto CLEANUP;

  for(int i = 0; i < clamp_edges.size(); i++){
    bounds[i] = 1.0; lower_or_upper[i] = 'U';
  }

  rval = PSEPlp_relaxbds(&m_lp, clamp_edges.size(), &clamp_edges[0],
			 &lower_or_upper[0], &bounds[0]);
  if(rval) goto CLEANUP;

  cout << "Relaxed bounds on: ";
  for(int i = 0; i < clamp_edges.size(); i++)
    cout << clamp_edges[i] << "\n";
  clamp_edges.clear();

 CLEANUP:
  if(rval){
    std::cerr << "Error entry point: AugHeuristic::clear_clamps()\n";
    clamp_edges.clear();
  }
  return rval;
}
