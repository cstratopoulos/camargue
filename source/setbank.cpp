#include "setbank.hpp"

#include<algorithm>
#include<iostream>
#include<iomanip>

using namespace std;

namespace CMR{

std::size_t hash_value(const IntervalSet &S){
  boost::hash<std::vector<IntPair>> hahser;
  return hahser(S.interval_list);
}

void interactive_test(){
  int node_count;
  cout << "Choose problem size: ";
  cin >> node_count;

  vector<int> tour_nodes(node_count);
  vector<int> perm(node_count);

  for(int i = 0; i < node_count; i++)
    tour_nodes[i] = i;

  random_shuffle(tour_nodes.begin(), tour_nodes.end());

  for(int i = 0; i < tour_nodes.size(); i++){
    perm[tour_nodes[i]] = i;
  }

  char answer = 'a';

  SetBank mybank(tour_nodes, perm);
  vector<HyperGraph> cut_bank;

  cout << "Testing ability to add and remove hypergraphs\n";
  while(answer == 'a' || answer == 'd'){
    if(answer == 'a'){
      int numsets;
      cout << "Enter number of sets in hypergraph: ";
      cin >> numsets;
      vector<vector<int>> node_sets(numsets);

      HyperGraph::CutType cut_t;

      if(numsets == 1)
	cut_t = HyperGraph::CutType::Segment;
      else if (numsets > 1)
	cut_t = HyperGraph::CutType::Blossom;
      else {
	cout << "Bad answer\n"; return;
      }
      
      for(int i = 0; i < numsets; i++){
	int setsize, current_node;
	cout << "Enter number of nodes in set " << i << ": ";
	cin >> setsize;
	cout << "Enter nodes in the set: ";
	for(int j = 0; j < setsize; j++){
	  cin >> current_node;
	  node_sets[i].push_back(current_node);
	}
      }

      cut_bank.emplace_back(HyperGraph(node_sets, cut_t));
    } else {
      int delete_ind;
      cout << cut_bank.size() << " cuts in bank. ";
      cout << "Which hypergraph do you want to delete? ";
      cin >> delete_ind;

      if(delete_ind < 0 || delete_ind >= cut_bank.size()){
	cout << "Index out of range\n";
      } else {
	(cut_bank.begin() + delete_ind)->delete_refs();
	cut_bank.erase(cut_bank.begin() + delete_ind);
      }
    }

    cout << "Keep testing? a to add a hypergraph, d to delete one, else quit: ";
    cin >> answer;
  }
}
}

using namespace CMR;

SetBank *HyperGraph::source_setbank;

IntervalSet::IntervalSet(vector<int> &nodelist,
			 const vector<int> &best_tour_nodes,
			 const vector<int> &perm) : use_count(0)
{
  sort(nodelist.begin(), nodelist.end(),
       [perm](const int node1, const int node2) -> bool {
	 return perm[node1] < perm[node2];
       });


  for(int i = 0; i < nodelist.size(); i++){
    int start_ind = perm[nodelist[i]];
    int end_ind = start_ind, cur_ind = start_ind;
    
    for(int j = i + 1; j < nodelist.size(); j++){
      cur_ind++;
      if(best_tour_nodes[cur_ind] == nodelist[j]){
	end_ind = cur_ind;
	i++;
      } else {
	break;
      }
    }
    
    interval_list.push_back(IntPair(start_ind, end_ind));
  }
}

SetBank::SetBank(vector<int> &best_tour_nodes,
		 vector<int> &_perm) :
  tour_nodes(best_tour_nodes),
  perm(_perm)
{
  HyperGraph::source_setbank = this;
}

int SetBank::extract_nodelist(const IntervalSet &interval_set,
			      vector<int> &result_nodelist){
  int rval = 0;
  result_nodelist.clear();
  const vector<IntPair> &interval_list = interval_set.interval_list;

  for(const IntPair &endpoints : interval_list){
    int start_ind = endpoints.first, end_ind = endpoints.second;

    try {
    for(int j = start_ind; j <= end_ind; j++)
      result_nodelist.push_back(tour_nodes[j]);
    } catch (const bad_alloc &){
      rval = 1; CMR_GOTO_CLEANUP("Couldn't push back result nodelist, ");
    }
  }

 CLEANUP:
  if(rval){
    cerr << "SetBank::extract_nodelist failed, result_nodelist invalid\n";
    result_nodelist.clear();
  }
  return rval;
}

IntervalSet *SetBank::add_or_increment(IntervalSet &newset)
{
  SetHash::iterator find_it = set_bank.find(newset);

  if(find_it != set_bank.end()){
    find_it->add_use();
  } else {
    newset.add_use();
    set_bank.insert(newset);
    find_it = set_bank.find(newset);
  }

  return const_cast<IntervalSet*>(&(*find_it));
}

void SetBank::del_or_decrement(IntervalSet &oldset)
{
  SetHash::iterator find_it = set_bank.find(oldset);

  if(find_it != set_bank.end()){
    find_it->del_use();
    if(find_it->use_count == 0){
      set_bank.erase(find_it);
    }
  }
}

HyperGraph::HyperGraph(vector<vector<int>> &node_sets,
		       const CutType _cut_type) :
  cut_type(_cut_type)
{
  for(vector<int> &current_set : node_sets){
    IntervalSet test_set(current_set, source_setbank->tour_nodes,
			 source_setbank->perm);
    IntervalSet *current_ref = source_setbank->add_or_increment(test_set);
    
    set_refs.push_back(current_ref);
  }
}

HyperGraph::HyperGraph(const IntPair &segment_ends) :
  cut_type(CutType::Segment)
{
  IntervalSet test_set{{segment_ends}};
  IntervalSet *new_ref = source_setbank->add_or_increment(test_set);
  set_refs.push_back(new_ref);
}

bool HyperGraph::same_tour(const vector<int> &best_tour_nodes)
{
  return best_tour_nodes == source_setbank->tour_nodes;
}

void HyperGraph::delete_refs(){
  for(IntervalSet *old_set : set_refs)
    source_setbank->del_or_decrement(*old_set);
}

