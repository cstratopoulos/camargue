#include "DPgraph.hpp"

#include <iostream>

using std::vector;
using std::cout;
using std::cerr;

namespace PSEP {

DPCutGraph::DPCutGraph(vector<vector<SimpleTooth::Ptr>> &_teeth,
		       const vector<int> &_perm,
		       const SupportGraph &_G_s) :
  light_teeth(_teeth), G_s(_G_s), perm(_perm) { CCcut_GHtreeinit(&gh_tree); }

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

  
  {
    int num_teeth = 0;
    for(const vector<SimpleTooth::Ptr> &vec : light_teeth)
      num_teeth += vec.size();
    
    if(cutgraph_nodes.size() != G_s.node_count + num_teeth + 1){
      cerr << "\tWRONG LIGHT TREE NODE COUNT!!!\n";
      rval = 1;
    }

    if(cut_ecap.size() != G_s.node_count + num_teeth){
      cerr << "\tWRONG LIGHT TREE EDGE COUNT!!!\n";
      rval = 1;
    }

    if(!rval)
      cout << "\tCORRECT TREE DIMS!!!!!\n"
	   << "\t ncount " << cutgraph_nodes.size() << ", ecount "
	   << cut_ecap.size() << "\n";
  }
  
 CLEANUP:
  if(rval)
    cerr << "DPCutGraph::build_light_tree failed\n";
  return rval;
}

int DPCutGraph::add_web_edges()
{
  int rval = 0;
  int start_ecount = cut_ecap.size();

  for(int end0 = 0; end0 < G_s.node_count; ++end0){
    for(int j = 0; j < G_s.nodelist[end0].s_degree; ++j){
      bool end0_in = false, end1_in = false;
      int end1 = G_s.nodelist[end0].adj_objs[j].other_end;
      
      if(end0 > end1) continue;

      double lp_weight = G_s.nodelist[end0].adj_objs[j].lp_weight;
      //search for smallest body with root end0 containing end1
      if(!light_teeth[end0].empty()){
	for(auto root_end0 = light_teeth[end0].begin(); //iterator
	    root_end0 != light_teeth[end0].end();
	    ++root_end0){
	  end1_in = (*root_end0)->body_contains(perm[end1]);
	  if(end1_in){
	    try { cut_elist.push_back((*root_end0)->cutgraph_index); }
	    catch (...) {
	      PSEP_SET_GOTO(rval, "Couldn't push back web edge. ");
	    }
	    break;
	  }
	}
      }
      if(!end1_in){ //none found, take degree eqn end0
	try { cut_elist.push_back(end0); } catch(...) {
	  PSEP_SET_GOTO(rval, "Couldn't push back web edge. ");
	}
      }

      if(!light_teeth[end1].empty()){
	for(auto root_end1 = light_teeth[end1].begin(); //iterator
	    root_end1 != light_teeth[end1].end();
	    ++root_end1){
	  end0_in = (*root_end1)->body_contains(perm[end0]);
	  if(end0_in){
	    try { cut_elist.push_back((*root_end1)->cutgraph_index); }
	    catch (...) {
	      PSEP_SET_GOTO(rval, "Couldn't push back web edge. ");
	    }
	    break;
	  }
	}
      }
      
      if(!end0_in){ //none found, take degree eqn end1
	try{ cut_elist.push_back(end1); } catch (...) {
	  PSEP_SET_GOTO(rval, "Couldn't push back web edge. ");
	}
      }

      try{ cut_ecap.push_back(lp_weight); } catch (...) {
	PSEP_SET_GOTO(rval, "Couldn't push back edge cap. ");
      }
    }
  }

  if(cut_ecap.size() != start_ecount + G_s.edge_count){
    rval = 1;
    cerr << "\t!!!WRONG WEB ECOUNT!!!!!\n";
  } else
    cout << "\t!!!CORRECT # WEB EDGES ADDED!!!!\n"
	 << "\tnew ecount: " << cut_ecap.size() << "\n";

  

 CLEANUP:
  if(rval)
    cerr << "DPCutGraph::add_web_edges failed\n";
  return rval;
}

int DPCutGraph::call_concorde_gomoryhu()
{
  int rval = 0;
  int ncount = cutgraph_nodes.size(), ecount = cut_ecap.size();
  int markcount = odd_nodes_list.size();

  int cutcount = 0;

  CCrandstate rstate;

  CCutil_sprand((int) real_zeit(), &rstate);

  cout << "\tCalling CCcut_gh with ncount: " << ncount << "\n"
       << "\tecount: " << ecount << "\n"
       << "\tmarkcount: " << markcount << ".......";
  rval = CCcut_gomory_hu(&gh_tree, ncount, ecount, &cut_elist[0], &cut_ecap[0],
			 markcount, &odd_nodes_list[0], &rstate);
  PSEP_CHECK_RVAL(rval, "CCcut_gomory_hu failed. ");
  cout << "Done.\n";

  cout << "\tPerforming odd cut dfs.....";
  dfs_odd_cuts(gh_tree.root, cutcount);
  cout << "Done, found " << cutcount << " odd cuts\n";

  

 CLEANUP:
  if(rval)
    cerr << "Problem in DPCutGraph::call_concorde_gomoryhu\n";
  return rval;
}

inline void DPCutGraph::dfs_odd_cuts(CC_GHnode *n, int &cutcount)
{
  if(n->parent)
    if(n->ndescendants % 2 == 1 && n->ndescendants > 1 && n->cutval < 0.9)
      ++cutcount;

  for(n = n->child; n; n = n->sibling)
    dfs_odd_cuts(n, cutcount);
}

}
