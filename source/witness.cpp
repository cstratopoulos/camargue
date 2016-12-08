#include "witness.hpp"
#include "err_util.hpp"

#include <stdexcept>
#include <string>

using std::vector;

using std::cout;
using std::cerr;
using std::string;

using std::exception;
using std::runtime_error;
using std::logic_error;

namespace CMR {


DPwitness::DPwitness(CandidateTeeth &cands,
                          const vector<int> &partition_nodes)
try :
  light_teeth(vector<vector<SimpleTooth>>(cands.light_teeth.size())),
  supp_dat(cands.supp_dat), perm(cands.best_dat.perm), CC_gh_q(1000)
  {
    CCcut_GHtreeinit(&gh_tree);

    for(int i : partition_nodes){
      int root = perm[i];
      for(const SimpleTooth::Ptr &T : cands.light_teeth[root])
        light_teeth[root].emplace_back(*T);
    }
  } catch(const exception &e) {
  cerr << e.what() << "\n";
  throw runtime_error("DPwitness constructor failed.");
 }

DPwitness::~DPwitness(){ CCcut_GHtreefree(&gh_tree); }

bool DPwitness::simple_DP_sep(CutQueue<dominoparity> &dp_q)
{
  try {
    build_light_tree();
    
    add_web_edges();
    
    build_gh_tree();
    if(CC_gh_q.empty()) return false;
    
    grab_dominos(dp_q);
  } catch(const exception &e){
    cerr << e.what() << "\n";
    throw runtime_error("Problem in DPwitness::simple_DP_sep.");
  }

  return (!dp_q.empty());
}

void DPwitness::build_light_tree()
{
  int num_cutnodes = 0, current_index, star_index;
  runtime_error err("Problem in DPwitness::build_light_tree. ");

  try { //nullptrs represent the degree eqn on each root
    for(int i = 0; i < light_teeth.size(); i++){
      cutgraph_nodes.push_back(nullptr);
    }
  } CMR_CATCH_PRINT_THROW("pushing back null nodes", err);

  current_index = cutgraph_nodes.size();

  try { //numbering and labelling cutgraph nodes
    for(vector<SimpleTooth> &t_vec : light_teeth)
      for(auto it = t_vec.rbegin(); it != t_vec.rend(); ++it){
	cutgraph_nodes.push_back(&(*it));
	it->cutgraph_index = current_index++;
      }
  } CMR_CATCH_PRINT_THROW("pushing back tooth nodes", err);


  cutgraph_nodes.push_back(nullptr);
  num_cutnodes = cutgraph_nodes.size();
  star_index = num_cutnodes - 1;
  //adding and labelling star node

  try { node_marks.resize(num_cutnodes, false); }
  CMR_CATCH_PRINT_THROW("initializing node marks bool", err);

  //adding degree eqn nodes between root and star
  try {//add edges between every root and central node for degree eqns
    for(int i = 0; i < light_teeth.size(); ++i){
      cut_elist.push_back(i);
      cut_elist.push_back(star_index);
      cut_ecap.push_back(0);
    }
  } CMR_CATCH_PRINT_THROW("Adding degree eqn nodes", err);

  try {
  for(int i = 0; i < light_teeth.size(); ++i){
    if(light_teeth[i].empty()) continue;

    for(auto child = light_teeth[i].begin();
	child != light_teeth[i].end(); ++child){
      int child_index = child->cutgraph_index;
      double child_slack = child->slack;
      bool found_parent = false;
      
      for(auto parent = child + 1; parent != light_teeth[i].end(); ++parent){
	found_parent = child->is_subset_of(*parent);
	if(found_parent){
	  int parent_index = parent->cutgraph_index;

	  cut_elist.push_back(child_index);
	  cut_elist.push_back(parent_index);
	  cut_ecap.push_back(child_slack);
	  
	  node_marks[child_index] = !(node_marks[child_index]);
	  node_marks[parent_index] = !(node_marks[parent_index]);

	  break;
	}
      }

      if(!found_parent){
	cut_elist.push_back(child_index);
	cut_elist.push_back(i);
	cut_ecap.push_back(child_slack);

	node_marks[child_index] = !(node_marks[child_index]);
	node_marks[i] = !(node_marks[i]);
      }
    }
  }
  } CMR_CATCH_PRINT_THROW("building light cut tree", err);

  try {
    for(int i = 0; i < node_marks.size(); i++)
      if(node_marks[i]){
	odd_nodes_list.push_back(i);
      }
  } CMR_CATCH_PRINT_THROW("setting gomoryhu marked nodes", err);

  
  int num_teeth = 0;
  for(const vector<SimpleTooth> &vec : light_teeth)
    num_teeth += vec.size();

  try { cg_delta_marks.resize(cutgraph_nodes.size(), 0); }
  CMR_CATCH_PRINT_THROW("allocating delta marks", err);  
}

void DPwitness::add_web_edges()
{
  runtime_error err("Problem in DPwitness::add_web_edges.");
  
  const vector<int> &support_elist = supp_dat.support_elist;
  const vector<double> &support_ecap = supp_dat.support_ecap;

  for(int i = 0; i < support_ecap.size(); ++i){
    double lp_weight = support_ecap[i];
    int end0 = perm[support_elist[2 * i]];
    int end1 = perm[support_elist[(2 * i) + 1]];
    int end0_container = -1, end1_container = -1;

    //search subtree root w  end0 for smallest (end0, S) with end1 in S
    for(const SimpleTooth &T : light_teeth[end0]){
      if(T.body_contains(end1)){
	end1_container = T.cutgraph_index;
	break;
      }
    }
    if(end1_container == -1) end1_container = end0;

    //search subtree w root end1 for smallest (end1, S) with end0 in S
    for(const SimpleTooth &T : light_teeth[end1]){
      if(T.body_contains(end0)){
	end0_container = T.cutgraph_index;
	break;
      }
    }
    if(end0_container == -1) end0_container = end1;

    try {
      cut_elist.push_back(end0_container);
      cut_elist.push_back(end1_container);
      cut_ecap.push_back(lp_weight);
    } CMR_CATCH_PRINT_THROW("pushing back web edge", err);
  }

  try {
    cutgraph_delta.resize(cut_ecap.size());
  } CMR_CATCH_PRINT_THROW("allocating cutgraph delta", err);
}

void DPwitness::build_gh_tree()
{
  runtime_error err("Problem in build_gh_tree.");
  int ncount = cutgraph_nodes.size();
  int ecount = cut_ecap.size();
  int markcount = odd_nodes_list.size();

  CCrandstate rstate;
  int seed = 99;
  
  CCutil_sprand(seed, &rstate);

  if(CCcut_gomory_hu(&gh_tree, ncount, ecount, &cut_elist[0], &cut_ecap[0],
                     markcount, &odd_nodes_list[0], &rstate)){
    throw runtime_error("CCcut_gomory_hu failed in DPwitness::build_gh_tree.");
  }

  try {
    dfs_odd_cuts(gh_tree.root);
  } CMR_CATCH_PRINT_THROW("dfsing tree for odd cuts", err);
}

void DPwitness::grab_dominos(CutQueue<dominoparity> &dp_q)
{
  runtime_error err("Problem in DPwitness::grab_dominos");
  int special_ind = cutgraph_nodes.size() - 1;

  while(!CC_gh_q.empty()){
    vector<int> cut_shore_nodes;
    int deltacount = 0;
    
    vector<int> handle_nodes;
    vector<SimpleTooth> used_teeth;
    vector<IntPair> nonneg_edges;

    try { expand_cut(CC_gh_q.peek_front(), cut_shore_nodes); }
    CMR_CATCH_PRINT_THROW("expanding cut nodes", err);

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
          CMR_CATCH_PRINT_THROW("adding nonneg edge", err);
	  continue;
	}

	try { used_teeth.push_back(*T1); }
        CMR_CATCH_PRINT_THROW("pushing back used tooth", err);

	continue;
      }

      if(!T1 && !T2){
	if(end0 != special_ind && end1 != special_ind){
	  try { nonneg_edges.push_back(IntPair(end0, end1)); }
          CMR_CATCH_PRINT_THROW("adding nonneg edge", err);

        continue;
      }

      try { handle_nodes.push_back(end0 == special_ind ? end1 : end0); }
      CMR_CATCH_PRINT_THROW("adding handle index", err);

      continue;
    }

      //at this point precisely one of T1, T2 is non-null
      try { used_teeth.push_back( T2 ? *T2 : *T1 ); }
      CMR_CATCH_PRINT_THROW("pushing back used tooth", err);
    }

    try {
      dominoparity newDP(used_teeth, handle_nodes, nonneg_edges);
      
      #pragma omp critical
      dp_q.push_back(newDP);
    } CMR_CATCH_PRINT_THROW("adding new dp ineq", err);

    CC_gh_q.pop_front();
  }
}

inline void DPwitness::dfs_odd_cuts(CC_GHnode *n)
{
  if(n->parent)
    if(n->ndescendants % 2 == 1 && n->ndescendants > 1 && n->cutval < 0.9) {
      if(CC_gh_q.empty() || n->cutval <= CC_gh_q.peek_front()->cutval)
	CC_gh_q.push_front(n);
    }

  for(n = n->child; n; n = n->sibling)
    dfs_odd_cuts(n);
}

inline void DPwitness::expand_cut(CC_GHnode *n, vector<int> &cut_nodes)
{
  for(int i = 0; i < n->listcount; ++i)
    cut_nodes.push_back(n->nlist[i]);

  for(n = n->child; n; n = n->sibling)
    expand_cut(n, cut_nodes);
}

}
