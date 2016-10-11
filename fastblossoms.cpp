#include "fastblossoms.hpp"
#include "PSEP_util.hpp"

#include <concorde/INCLUDE/cut.h>

#include <iostream>

using namespace std;
using namespace PSEP;

int Cut<fastblossom>::separate()
{
  int rval = 1;
  int ncount = m_graph.node_count;
  int component_count = 0;

  //array of length component_count containing sizes of components
  int *component_sizes = (int *) NULL;

  //array of length ncount storing component nodelists consecutively
  //the first component_sizes[0] elements making up the first component, etc
  int *component_nodes = (int *) NULL;

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
    } else {
      //cout << "Not a fast blossom : (" << endl;
    }
  }
  }

  if(local_q.empty()) rval = 2;


 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<fastblossom>::separate\n";
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
  frac_indices.clear();
  frac_elist.clear();
  frac_ecap.clear();
  return rval;
}
