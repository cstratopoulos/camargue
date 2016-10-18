#include "fastblossoms.hpp"
#include "PSEP_util.hpp"

#include <concorde/INCLUDE/cut.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <iterator>

#include <cmath>

using namespace std;
using namespace PSEP;

int Cut<fastblossom>::oc_sep()
{
  int rval = 0;
  int ncount = m_graph.node_count;
  int component_count = 0;

  //array of length component_count containing sizes of components
  int *component_sizes = (int *) NULL;

  //array of length ncount storing component nodelists consecutively
  //the first component_sizes[0] elements making up the first component, etc
  int *component_nodes = (int *) NULL;

  vector<int> frac_indices;
  vector<int> frac_elist;
  vector<double> frac_ecap;

  try {
    for(int i = 0; i < support_indices.size(); i++){
      int edge_index = support_indices[i];
      if(support_ecap[i] <= 1 - LP::EPSILON){
	// cout << "Adding edge " << edge_index << ", ends "
	//      << support_elist[2*i] << ", "
	//      << support_elist[2*i + 1] << endl;
	frac_indices.push_back(edge_index);
	frac_ecap.push_back(support_ecap[i]);
	frac_elist.push_back(support_elist[2 * i]);
	frac_elist.push_back(support_elist[(2 * i) + 1]);
      }
    }
  } catch (...) { rval = 1; PSEP_GOTO_CLEANUP("Couldn't set frac edges, "); }

  if(frac_indices.empty()){
    //    cout << "Frac graph is empty!" << endl;
    rval = 2; goto CLEANUP;
  }

  rval = CCcut_connect_components(ncount, frac_indices.size(),
				  &frac_elist[0], &frac_ecap[0],
				  &component_count,
				  &component_sizes,
				  &component_nodes);
  PSEP_CHECK_RVAL(rval, "CCcut_connect_components failed, ");

  
  {int j = 0; //scoped index for traversing component nodes
  for(int i = 0; i < component_count; i++){
    if(component_sizes[i] == 1){
      j++;
      continue;
    }
    //cout << "Component " << i << " has size " << component_sizes[i] << endl;
    vector<int> handle_nodes, edge_indices;
    int deltacount = 0;
    bool has_intersection = false;

    try{
      for(int l = 0; l < component_sizes[i]; l++){
	handle_nodes.push_back(component_nodes[j]);
	//	cout << component_nodes[j] << endl;
	j++;
      }
    } catch(...) { rval = 1; PSEP_GOTO_CLEANUP("Couldn't get handle nodes, ");}

    // cout << "Stored handle in vector" << endl;

    do {
      GraphUtils::get_delta(handle_nodes.size(), &handle_nodes[0],
			    support_indices.size(), &support_elist[0],
			    &deltacount, &delta[0], &edge_marks[0]);
      //cout << "Got delta, deltacount: " << deltacount << ", edges:" << endl;
      edge_indices.clear();
      try {
	for(int k = 0; k < deltacount; k++){
	  int sup_ind = delta[k];
	  int edge_index = support_indices[sup_ind];
	  if(support_ecap[sup_ind] > 1 - LP::EPSILON){
	    edge_marks[support_elist[2*sup_ind]] += 1;
	    edge_marks[support_elist[(2*sup_ind) + 1]] += 1;
	    //cout << "Pushing back " << edge_index << endl;
	    edge_indices.push_back(edge_index);
	  }
	}
      } catch(...){ rval = 1; PSEP_GOTO_CLEANUP("Couldn't grow teeth. "); }

      // cout << "Got collection of " << edge_indices.size() << " candiates"
      // 	   << endl;
      // cout << "Testing for intersection..." << endl;
      
      has_intersection = false;
      try {
	for(int l = 0; l < edge_marks.size(); l++)
	  if(edge_marks[l] > 1){
	    // cout << "Node " << l << " meets " << edge_marks[l] << ", "
	    // 	 << "adding it to handle" << endl;
	    handle_nodes.push_back(l);
	    has_intersection = true;
	  }
      } catch(...) {
	rval = 1; PSEP_GOTO_CLEANUP("Couldn't grow handle nodes. ");
      }

      for(int k = 0; k < deltacount; k++){
	int sup_ind = delta[k];
	if(support_ecap[sup_ind] > 1 - LP::EPSILON){
	  edge_marks[support_elist[2*sup_ind]] = 0;
	  edge_marks[support_elist[(2*sup_ind) + 1]] = 0;
	}
      }
    } while(edge_indices.size() > 3 && has_intersection);

    // cout << "    Out of do-while loop. Intersection: " << has_intersection
    // 	 << ", num teeth: " << edge_indices.size() << endl;

    if(!has_intersection && edge_indices.size() >= 3 &&
       edge_indices.size() % 2 == 1){
      // cout << "    ****Found fast blossom!**** "
      // 	   << "Testing if it is tight at best tour..." << endl;
      bool tight = false;
      int num_teeth = edge_indices.size();
      //      cout << "|F| = " << num_teeth << ", ";
      
      int sum_e_F = 0;
      int sum_handle_minus_F = 0;

      for(int k = 0; k < num_teeth; k++){
	sum_e_F += best_tour_edges[edge_indices[k]];
	sum_handle_minus_F -= best_tour_edges[edge_indices[k]];
      }

      //cout << "x(F) = " << sum_e_F << ", ";
      
      GraphUtils::get_delta(handle_nodes, m_graph.edges, &deltacount,
			    delta, edge_marks);

      for(int k = 0; k < deltacount; k++){
	sum_handle_minus_F += best_tour_edges[delta[k]];
      }

      //cout << "x(d(H)\\F) = " << sum_handle_minus_F << endl;

      tight = ((sum_e_F == num_teeth && sum_handle_minus_F == 1) ||
	       (sum_e_F == num_teeth - 1 && sum_handle_minus_F == 0));

      // cout << "    Cut is " << (tight ? "TIGHT!!!" : "not tight :(")
      // 	   << endl;
      if(tight){
	fastblossom newcut(handle_nodes, edge_indices);
	local_q.push_front(newcut);
	//cout << "    Pushed blossom to local queue" << endl;
      }
    }
  }
  }

  if(local_q.empty()) rval = 2;


 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<fastblossom>::oc_sep\n";
  if(component_nodes) free(component_nodes);
  if(component_sizes) free(component_sizes);
  return rval;
}

int Cut<fastblossom>::GH_sep()
{
  cout << "Calling GH sep-------" << endl;
  int rval = 0;
  int ncount = m_graph.node_count;
  int component_count = 0;

  //array of length component_count containing sizes of components
  int *component_sizes = (int *) NULL;

  //array of length ncount storing component nodelists consecutively
  //the first component_sizes[0] elements making up the first component, etc
  int *component_nodes = (int *) NULL;

  vector<int> eps_indices;
  vector<int> eps_elist;
  vector<double> eps_ecap;

  try {
    for(int i = 0; i < support_indices.size(); i++){
      int edge_index = support_indices[i];
      if(support_ecap[i] >= GH_eps && support_ecap[i] <= 1 - GH_eps){
	eps_indices.push_back(edge_index);
	eps_ecap.push_back(support_ecap[i]);
	eps_elist.push_back(support_elist[2 * i]);
	eps_elist.push_back(support_elist[(2 * i) + 1]);
      }
    }
  } catch(...){ rval = 1; PSEP_GOTO_CLEANUP("Couldn't set eps edges. "); }

  if(eps_indices.empty()) {
    rval = 2;
    cout << "No eps edges!" << endl;
    goto CLEANUP;
  }

  cout << "Got nonempty collection of epsilon edges" << endl;

  rval = CCcut_connect_components(ncount, eps_indices.size(),
				  &eps_elist[0], &eps_ecap[0],
				  &component_count,
				  &component_sizes,
				  &component_nodes);
  PSEP_CHECK_RVAL(rval, "CCcut_connect_components failed. ");

  cout << "Got components of gh cutgraph" << endl;

  {int j = 0; //scoped index for traversing component nodes
    for(int i = 0; i < component_count; i++){
      if(component_sizes[i] == 1){
	j++;
	continue;
      }

      cout << "Nontrivial component of size" << component_sizes[i] << endl;

      vector<int> handle_nodes;
      vector<int> teeth;
      int deltacount = 0;
      bool has_intersection = true;

      try {
	for(int k = 0; k < component_sizes[i]; k++){
	  handle_nodes.push_back(component_nodes[j]);
	  j++;
	}
      } catch(...){
	rval = 1; PSEP_GOTO_CLEANUP("Couldn't get handle nodes. ");
      }

      cout << "Got " << handle_nodes.size() << " handle nodes" << endl;

      do {
	cout << "Pass of do while loop" << endl;
	teeth.clear();
	GraphUtils::get_delta(handle_nodes.size(), &handle_nodes[0],
			      support_indices.size(), &support_elist[0],
			      &deltacount, &delta[0], &edge_marks[0]);
	cout << "Got delta" << endl;

	try {
	  for(int k = 0; k < deltacount; k++){
	    int sup_ind = delta[k];
	    int edge_ind = support_indices[sup_ind];

	    if(support_ecap[sup_ind] > 1 - GH_eps){
	      edge_marks[support_elist[2 * sup_ind]] += 1;
	      edge_marks[support_elist[(2 * sup_ind) + 1]] += 1;
	      teeth.push_back(edge_ind);
	    }
	  }
	} catch(...){ rval = 1; PSEP_GOTO_CLEANUP("Couldn't get teeth. "); }

	cout << "Got " << teeth.size() << " teeth" << endl;

	if(teeth.size() % 2 == 0){
	  int best_edge_ind = -1, best_sup_ind = -1;
	  double max_little = -1.0;

	  for(int k = 0; k < deltacount; k++){
	    int sup_ind = delta[k];
	    int edge_ind = support_indices[sup_ind];

	    if(support_ecap[sup_ind] < GH_eps &&
	       support_ecap[sup_ind] > max_little){
	      best_edge_ind = edge_ind;
	      best_sup_ind = sup_ind;
	      max_little = support_ecap[sup_ind];
	    }
	  }

	  if(best_edge_ind >= 0){
	    edge_marks[support_elist[2 * best_sup_ind]] += 1;
	    edge_marks[support_elist[(2 * best_sup_ind) + 1]] += 1;
	    teeth.push_back(best_edge_ind);
	    cout << "Now have " << teeth.size() << " teeth" << endl;
	  } else {
	    break;
	    cout << "Couldn't fix teeth size!" << endl;
	  }
	}

	has_intersection = false;

	for(int k = 0; k < handle_nodes.size(); k++)
	  edge_marks[handle_nodes[k]] *= -1;

	for(int k = 0; k < edge_marks.size(); k++){
	  if(edge_marks[k] > 1){
	    try { handle_nodes.push_back(k); } catch(...){
	      rval = 1; PSEP_GOTO_CLEANUP("Couldn't grow handle. ");
	    }

	    cout << "Added node " << k << ", meets out of handle" << endl;
	    has_intersection = true;
	    continue;
	  }

	  if(edge_marks[k] < -1){
	    handle_nodes.erase(remove(handle_nodes.begin(),
				      handle_nodes.end(), k),
			       handle_nodes.end());
	    cout << "Removed " << k << ", intersection in handle" << endl;
	    has_intersection = true;
	    continue;
	  }	  
	}

	for(int k = 0; k < edge_marks.size(); k++)
	  edge_marks[k] = 0;
      
      } while(has_intersection);

      cout << "Out of do while loop, checking for even number of teeth or "
	   << "less than 3" << endl;
      if(teeth.size() % 2 == 0 || teeth.size() < 3){
	cout << "Bad teeth size, continuing" << endl;
	continue;
      }
      
      
      sort(teeth.begin(), teeth.end(),
	   [this](const int &e1, const int &e2) -> bool
	   { return m_lp_edges[e1] > m_lp_edges[e2]; });

      cout << "Sorted teeth by lpweight" << endl;

      int k_upper = teeth.size();
      int k_lower = (k_upper == 3) ? k_upper : k_upper - 2;

      GraphUtils::get_delta(handle_nodes, m_graph.edges, &deltacount, delta,
			    edge_marks);

      cout << "Got handle delta, k lower/upper: "
	   << k_lower << ", " << k_upper << endl;

      for(int k = k_lower; k <= k_upper; k += 2){
	bool tight = false, violated = false;
	int num_teeth = k;
	cout << "Num teeth (k) = " << k << endl;

	int sum_e_F = 0, sum_handle_minus_F = 0;

	double lp_e_F = 0.0, lp_handle_minus_F = 0.0;
	double rhs = 1 - num_teeth, lhs;

	for(int l = 0; l < num_teeth; l++){
	  sum_e_F += best_tour_edges[teeth[l]];
	  lp_e_F += m_lp_edges[teeth[l]];

	  sum_handle_minus_F -= best_tour_edges[teeth[l]];
	  lp_handle_minus_F -= m_lp_edges[teeth[l]];
	}

	cout << "x(F) tour/lp: " << sum_e_F << ", " << lp_e_F << endl;

	for(int l = 0; l < deltacount; l++){
	  sum_handle_minus_F += best_tour_edges[delta[l]];
	  lp_handle_minus_F += m_lp_edges[delta[l]];
	}

	cout << "x(d(H)\\F) tour/lp: " << sum_handle_minus_F
	     << ", " << lp_handle_minus_F << endl;

	tight = ((sum_e_F == num_teeth && sum_handle_minus_F == 1) ||
		 (sum_e_F == num_teeth - 1 && sum_handle_minus_F == 0));

	lhs = lp_handle_minus_F - lp_e_F;
	violated = (lhs < rhs);

	cout << "tight/viol " << tight << "/" << violated << endl;

	if(!(tight && violated)) continue;

	try {
	  if(num_teeth == teeth.size()){
	    local_q.push_front(fastblossom(handle_nodes, teeth));
	  } else {
	    vector<int> partial_teeth(teeth.begin(),
				      teeth.begin() + num_teeth);
	    local_q.push_front(fastblossom(handle_nodes, partial_teeth));
	  }
	} catch(...){
	  rval = 1; PSEP_GOTO_CLEANUP("Couldn't push to local queue. ");
	}
	cout << "Pushed to local q" << endl;
      }      
    }
    }

  if(local_q.size() == 0) rval = 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Cut<fastblossom>::GH_sep failed.\n";
  if(component_nodes) free(component_nodes);
  if(component_sizes) free(component_sizes);

  return rval;
}

int Cut<fastblossom>::build_hypergraph(const fastblossom &blossom_cut)
{
  int rval = 0;
  vector<vector<int>> node_sets;

  try { node_sets.push_back(blossom_cut.handle); } catch(...){
    rval = 1; PSEP_GOTO_CLEANUP("Couldn't push back handle. ");
  }

  try {
    for(int i = 0; i < blossom_cut.edge_indices.size(); i++){
      int current_index = blossom_cut.edge_indices[i];
      vector<int> new_tooth{m_graph.edges[current_index].end[0],
			    m_graph.edges[current_index].end[1]};
      node_sets.push_back(new_tooth);
    }
  } catch(...){ rval = 1; PSEP_GOTO_CLEANUP("Couldn't push back teeth. "); }

  try {
    HyperGraph newcut(node_sets, HyperGraph::CutType::Blossom);
    external_q.push_back(newcut);
  } catch(...){ rval = 1; PSEP_GOTO_CLEANUP("Couldn't push hypergraph. "); }

 CLEANUP:
  if(rval)
    cerr << "Cut<fastblossom>::build_hypergraph failed\n";
  return rval;
}

int Cut<fastblossom>::separate()
{
  int rval = 0, oc_rval = 0, gh_rval = 2;

  oc_rval = oc_sep();
  if(oc_rval == 1) { rval = 1; goto CLEANUP; }

  if(oc_rval == 2){
    cout << "No odd cut blossoms found, calling GH" << endl;
    gh_rval = GH_sep();
    if(gh_rval == 1) { rval = 1; goto CLEANUP; }
  }

  if(gh_rval == 2 && oc_rval == 2) rval = 2;

  if(gh_rval == 1) rval = 1;

 CLEANUP:
  if(rval == 1)
    cerr << "Cut<fastblossom>::separate failed\n";
  return rval;
}

//TODO: this should be fully templatized
int Cut<fastblossom>::add_cuts()
{
  int rval = 0;

  while(!local_q.empty()){
    rval = build_hypergraph(local_q.peek_front());
    if(rval) goto CLEANUP;

    local_q.pop_front();
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in Cut<fastblossom>::add_cuts\n";
  return rval;
}

int Cut<fastblossom>::cutcall()
{
  //  cout << "--------------------The fastblossom cutcall" << endl;
  int rval = separate();
  if(rval) goto CLEANUP;

  rval = add_cuts();
  if(rval) goto CLEANUP;
  // cout << "Added cuts to external q, now has size "
  //      << external_q.size() << endl;
  
 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cuts<fastblossom>::cutcall\n";
  return rval;
}
