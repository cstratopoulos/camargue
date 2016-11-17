#include "DPgraph.hpp"
#include "PSEP_util.hpp"
#include "Graph.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>

using std::vector;
using std::cout;
using std::cerr;

#define PSEP_DO_VIZ

namespace PSEP {

DPCutGraph::DPCutGraph(vector<vector<SimpleTooth::Ptr>> &_teeth,
		       CandidateTeeth &_cands,
		       const vector<int> &_perm,
		       const SupportGraph &_G_s) :
  light_teeth(_teeth),
  cands(_cands),
  G_s(_G_s),
  perm(_perm),
  CC_gh_q(25)
{
  CCcut_GHtreeinit(&gh_tree);
}

DPCutGraph::~DPCutGraph(){ CCcut_GHtreefree(&gh_tree); }

int DPCutGraph::simple_DP_sep(CutQueue<dominoparity> &domino_q)
{
  int rval = 0;

  rval = build_light_tree();
  if(rval) goto CLEANUP;

  rval = add_web_edges();
  if(rval) goto CLEANUP;

  rval = call_concorde_gomoryhu();
  if(rval) goto CLEANUP;

  rval = grab_cuts(domino_q);
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in DPCutGraph::simple_DP_sep.\n";
  return rval;
}

int DPCutGraph::build_light_tree()
{
  int rval = 0;
  int num_cutnodes = 0, current_index, star_index;

#ifdef PSEP_DO_VIZ
  std::ofstream light_out;
  light_out.open(ofname + "-light.gv");
  light_out << "graph {\n\t"
	    << "mindist=0.5;\n\t"; //circo
#endif


  try { //nullptrs represent the degree eqn on each root
    for(int i = 0; i < light_teeth.size(); i++){
      cutgraph_nodes.push_back(nullptr);
#ifdef PSEP_DO_VIZ
      light_out << i << "[color=aquamarine3];\n\t";
#endif
    }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't push back degree eqn nodes. "); }

  current_index = cutgraph_nodes.size();

#ifdef PSEP_DO_VIZ
  light_out << "//numbering and labelling cutgraph node for each tooth\n\t";
#endif
  try {
    for(const vector<SimpleTooth::Ptr> &t_vec : light_teeth)
      for(auto it = t_vec.rbegin(); it != t_vec.rend(); ++it){
	cutgraph_nodes.push_back(it->get());
	(*it)->cutgraph_index = current_index++;
#ifdef PSEP_DO_VIZ
	light_out << (*it)->cutgraph_index
		  << "[label=\"" << cands.print_label(**it, true)
		  << "\"];\n\t";
#endif
      }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't add tooth ineqs. " ); }


  cutgraph_nodes.push_back(nullptr);
  num_cutnodes = cutgraph_nodes.size();
  star_index = num_cutnodes - 1;
#ifdef PSEP_DO_VIZ
  light_out << "//adding and labelling star node\n\t";
  light_out << star_index << "[label=\"X\"];\n\t";
  light_out << "root=" << star_index << ";\n\t";
#endif

  try { node_marks.resize(num_cutnodes, false); } catch(...){
    PSEP_SET_GOTO(rval, "Couldn't initialize cut marks bool. ");
  }

#ifdef PSEP_DO_VIZ
  light_out << "//adding degree eqn nodes between root and star\n\t";
#endif
  try {//add edges between every root and central node for degree eqns
    for(int i = 0; i < light_teeth.size(); ++i){
      cut_elist.push_back(i);
      cut_elist.push_back(star_index);
      cut_ecap.push_back(0);
#ifdef PSEP_DO_VIZ
      light_out << star_index << "--" << i << "[label=\"d" << i << "\"];\n\t";
#endif
    }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't add degree eqn edges. "); }

  try {
  for(int i = 0; i < light_teeth.size(); ++i){
    if(light_teeth[i].empty()) continue;

    for(auto child = light_teeth[i].begin();
	child != light_teeth[i].end(); ++child){
      int child_index = (*child)->cutgraph_index;
      double child_slack = (*child)->slack;
      bool found_parent = false;
      
      for(auto parent = child + 1; parent != light_teeth[i].end(); ++parent){
	found_parent = (*child)->is_subset_of(**parent);
	if(found_parent){
	  int parent_index = (*parent)->cutgraph_index;

	  cut_elist.push_back(child_index);
	  cut_elist.push_back(parent_index);
#ifdef PSEP_DO_VIZ
	  light_out << child_index << "--" << parent_index
		    << "[label=\"" << std::setprecision(2)
		    << child_slack << "\"]"
		    << ";\n\t";
#endif	  
	  cut_ecap.push_back(child_slack);
	  
	  node_marks[child_index] = !(node_marks[child_index]);
	  node_marks[parent_index] = !(node_marks[parent_index]);

	  break;
	}
      }

      if(!found_parent){
	cut_elist.push_back(child_index);
	cut_elist.push_back(i);
#ifdef PSEP_DO_VIZ
	light_out << child_index << "--" << i
		  << "[label=\"" << std::setprecision(2)
		  << child_slack << "\"]"
		  << ";\n\t";
#endif
	cut_ecap.push_back(child_slack);

	node_marks[child_index] = !(node_marks[child_index]);
	node_marks[i] = !(node_marks[i]);
      }
    }
  }
  } catch(...) { PSEP_SET_GOTO(rval, "Problem building light cuttree. "); }

  try {
    for(int i = 0; i < node_marks.size(); i++)
      if(node_marks[i]){
	odd_nodes_list.push_back(i);
#ifdef PSEP_DO_VIZ
	light_out << i << "[color=red];\n\t";
#endif
      }
  } catch(...){ PSEP_SET_GOTO(rval, "Couldn't set gomoryhu marked nodes. "); }

  
  {
    int num_teeth = 0;
    for(const vector<SimpleTooth::Ptr> &vec : light_teeth)
      num_teeth += vec.size();
    
    if(cutgraph_nodes.size() != G_s.node_count + num_teeth + 1){
      PSEP_SET_GOTO(rval, "Wrong light tree node count. ");
    }

    if(cut_ecap.size() != G_s.node_count + num_teeth){
      PSEP_SET_GOTO(rval, "Wrong light tree edge count. ");
    }
  }

  try { cg_delta_marks.resize(cutgraph_nodes.size(), 0); } catch (...) {
    PSEP_SET_GOTO(rval, "Couldn't allocate delta node marks. ");
  }
  
 CLEANUP:
  if(rval)
    cerr << "DPCutGraph::build_light_tree failed\n";
#ifdef PSEP_DO_VIZ
  light_out << "}";
  light_out.close();
#endif
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
    PSEP_SET_GOTO(rval, "Wrong number of web edges added. ");
  } 

  try { cutgraph_delta.resize(cut_ecap.size()); } catch (...) {
    PSEP_SET_GOTO(rval, "Couldn't allocate cutgraph delta. ");
  }
  

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

  CCrandstate rstate;
  //int seed = (int) real_zeit();
  CCutil_sprand(1479252604, &rstate);

  cout << "\tCalling CCcut_gh with ncount: " << ncount << "\n"
       << "\tecount: " << ecount << "\n"
       << "\tmarkcount: " << markcount << ".......";
  rval = CCcut_gomory_hu(&gh_tree, ncount, ecount, &cut_elist[0], &cut_ecap[0],
			 markcount, &odd_nodes_list[0], &rstate);
  PSEP_CHECK_RVAL(rval, "CCcut_gomory_hu failed. ");
  cout << "Done.\n";

  cout << "\tPerforming odd cut dfs.....";
  
  try { dfs_odd_cuts(gh_tree.root); } catch (...) {
    PSEP_SET_GOTO(rval, "Problem dfsing tree for odd cuts. ");
  }
  
  cout << "Done. " << CC_gh_q.size() << " cuts in queue. ";
  if(CC_gh_q.size())
    cout << "Best cutval: " << CC_gh_q.peek_front()->cutval;
  cout << "\n";

  if(CC_gh_q.empty()) rval = 2;

  CCcut_GHtreeprint(&gh_tree);

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in DPCutGraph::call_concorde_gomoryhu\n";
  return rval;
}

inline void DPCutGraph::dfs_odd_cuts(CC_GHnode *n)
{
  if(n->parent)
    if(n->ndescendants % 2 == 1 && n->ndescendants > 1 && n->cutval < 0.9) {
      if(CC_gh_q.empty() || n->cutval < CC_gh_q.peek_front()->cutval)
	CC_gh_q.push_front(n);
      else
	CC_gh_q.push_back(n);
    }

  for(n = n->child; n; n = n->sibling)
    dfs_odd_cuts(n);
}

inline void DPCutGraph::expand_cut(CC_GHnode *n, vector<int> &cut_nodes)
{
  for(int i = 0; i < n->listcount; ++i)
    cut_nodes.push_back(n->nlist[i]);

  for(n = n->child; n; n = n->sibling)
    expand_cut(n, cut_nodes);
}

int DPCutGraph::grab_cuts(CutQueue<dominoparity> &domino_q)
{
  int rval = 0;
  int special_ind = cutgraph_nodes.size() - 1;

  while(!CC_gh_q.empty()){
    vector<int> cut_shore_nodes;
    int deltacount = 0;
    
    vector<int> handle_nodes;
    vector<SimpleTooth*> used_teeth;
    vector<IntPair> nonneg_edges;

    try { expand_cut(CC_gh_q.peek_front(), cut_shore_nodes); } catch (...) {
      PSEP_SET_GOTO(rval, "Couldn't expand cut nodes. ");
    }

    GraphUtils::get_delta(cut_shore_nodes.size(), &cut_shore_nodes[0],
			  cut_ecap.size(), &cut_elist[0], &deltacount,
			  &cutgraph_delta[0], &cg_delta_marks[0]);

    for(int i = 0; i < deltacount; ++i){
      int edge_ind = cutgraph_delta[i];
      int end0 = cut_elist[2 * edge_ind], end1 = cut_elist[(2 * edge_ind) + 1];

      SimpleTooth *T1 = cutgraph_nodes[end0];
      SimpleTooth *T2 = cutgraph_nodes[end1];

      if(T1 && T2){
	if(T1->root != T2->root){
	  try { nonneg_edges.push_back(IntPair(T1->root, T2->root)); }
	  catch (...) { PSEP_SET_GOTO(rval, "Couldn't add nonneg edge. "); }
	  continue;
	}

	try { used_teeth.push_back(T1); } catch (...) {
	  PSEP_SET_GOTO(rval, "Couldn't push back used tooth. ");
	}

	continue;
      }

      if(!T1 && !T2){
	if(end0 != special_ind && end1 != special_ind){
	  try { nonneg_edges.push_back(IntPair(end0, end1)); } catch (...) {
	    PSEP_SET_GOTO(rval, "Couldn't add nonneg edge. ");
	  }

	  continue;
	}

	try { handle_nodes.push_back(end0 == special_ind ? end1 : end0); }
	catch (...) { PSEP_SET_GOTO(rval, "Couldn't add handle ind. "); }

	continue;
      }

      //at this point precisely one of T1, T2 is non-null
      try { used_teeth.push_back( T2 ? T2 : T1 ); } catch (...) {
	PSEP_SET_GOTO(rval, "Couldn't push back used tooth. ");
      }
    }

    try {
      dominoparity newDP(used_teeth, handle_nodes, nonneg_edges);
      domino_q.push_back(newDP);
    } catch (...) { PSEP_SET_GOTO(rval, "Couldn't push back new DP ineq. "); }

    CC_gh_q.pop_front();
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in DPCutGraph::grab_cuts.\n";
  return rval;
}



}
