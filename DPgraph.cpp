#include "DPgraph.hpp"

#include <iostream>

using std::vector;

namespace PSEP {

DPCutGraph::DPCutGraph(const vector<vector<SimpleTooth::Ptr>> &_teeth) :
  light_teeth(_teeth) { CCcut_GHtreeinit(&gh_tree); }

DPCutGraph::~DPCutGraph(){ CCcut_GHtreefree(&gh_tree); }

int DPCutGraph::build_light_tree()
{
  int rval = 0;
  int num_cutnodes = 0, current_index, star_index;


  try { //nullptrs represent the degree eqn on each root
    for(int i = 0; i < light_teeth.size(); i++)
      cutgraph_nodes.push_back(nullptr);
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't push back degree eqn nodes. "); }

  current_index = cutgraph_nodes.size();

  try { //number a cutgraph node for each tooth ineq, add to vec
    for(const vector<SimpleTooth::Ptr> &t_vec : light_teeth)
      for(auto it = t_vec.rbegin(); it != t_vec.rend(); ++it){
	cutgraph_nodes.push_back(it->get());
	(*it)->cutgraph_index = current_index++;
      }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't add tooth ineqs. " ); }

  //add a node representing the special central node, "star" node
  cutgraph_nodes.push_back(nullptr);
  num_cutnodes = cutgraph_nodes.size();
  star_index = num_cutnodes - 1;

  try { node_marks.resize(num_cutnodes, false); } catch(...){
    PSEP_SET_GOTO(rval, "Couldn't initialize cut marks bool. ");
  }

  try {//add edges between every root and central node for degree eqns
    for(int i = 0; i < light_teeth.size(); ++i){
      cut_elist.push_back(i);
      cut_elist.push_back(star_index);
      cut_ecap.push_back(0);
    }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't add degree eqn edges. "); }

  try {
    for(int i = 0; i < light_teeth.size(); i++){
      if(light_teeth[i].empty()) continue;

      auto child = light_teeth[i].rbegin();
      int child_index = (*child)->cutgraph_index;
      double child_slack = (*child)->slack;

      cut_elist.push_back(child_index);
      cut_elist.push_back(i);
      cut_ecap.push_back(child_slack);

      //complement bits for odd or evenness
      node_marks[i] = !(node_marks[i]);
      node_marks[child_index] = !(node_marks[child_index]);

      ++child; //for loop initialized one after begin
      for(; child != light_teeth[i].rend(); ++child){
	auto parent = child;
	int parent_index = (*parent)->cutgraph_index;
	bool found_parent;
	
	while(true){
	  --parent;

	  found_parent = (**child).is_subset_of(**parent);

	  if(found_parent){//add edge between smallest parent
	    cut_elist.push_back(child_index);
	    cut_elist.push_back(parent_index);
	    cut_ecap.push_back(child_slack);

	    node_marks[child_index] = !(node_marks[child_index]);
	    node_marks[parent_index] = !(node_marks[parent_index]);

	    break;
	  }

	  if(parent == light_teeth[i].rbegin()) break;
	}

	if(!found_parent){//then edge is between degree node
	  cut_elist.push_back(child_index);
	  cut_elist.push_back(i);
	  cut_ecap.push_back(child_slack);

	  node_marks[i] = !(node_marks[i]);
	  node_marks[child_index] = !(node_marks[child_index]);
	}	
      }      
    }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't add tree edges. "); }

  try{
    for(int i = 0; i < node_marks.size(); i++)
      if(node_marks[i])
	odd_nodes_list.push_back(i);
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't set gomoryhu marked nodes. "); }

 CLEANUP:
  if(rval)
    std::cerr << "DPCutGraph::build_light_tree failed\n";
  return rval;
}

int DPCutGraph::add_web_edges()
{
  int rval = 0;

  //CLEANUP:
  if(rval)
    std::cerr << "DPCutGraph::add_web_edges failed\n";
  return rval;
}

}
