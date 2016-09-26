#include "setbank.hpp"

#include<algorithm>
#include<iostream>
#include<iomanip>

using namespace std;

namespace PSEP{

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

  SetHash set_bank;
  cout << "Declared empty set_bank\n";

  char answer = 'y';

  while(answer == 'y'){
    vector<int> set_nodes;
    int set_size, set_elem;
    cout << "Enter size of nodeset to test: ";
    cin >> set_size;
    cout << "Enter nodes in set: \n";
    for(int i = 0; i < set_size; i++){
      cin >> set_elem;
      set_nodes.push_back(set_elem);
    }

    IntervalSet newset(set_nodes, tour_nodes, perm);
    cout << "Address of newset: " << &newset << "\n";

    SetHash::iterator get = set_bank.find(newset);
    if(get == set_bank.end()){
      cout << "Adding new element, address in set: ";
      newset.add_use();
      set_bank.insert(newset);
      get = set_bank.find(newset);
      cout << &(*get) << "\n";
    } else {
      cout << "Incrementing existing element, address: "
	   << &(*get) << "\n";
      get->add_use();
      
    }
    
    cout << "Size of set bank is now: " << set_bank.size() << "\n";

    cout << "Continue? y for yes: ";
    cin >> answer;
  }
  cout << "Done testing, printing use counts/addresses:\n";
  for(SetHash::iterator it = set_bank.begin(); it != set_bank.end(); it++){
    cout << it->use_count << ", " << &(*it) << "\n";
  }
}

}

using namespace PSEP;

unique_ptr<SetBank> HyperGraph::source_setbank;

IntervalSet::IntervalSet(vector<int> &nodelist,
			 const vector<int> &best_tour_nodes,
			 const vector<int> &perm) : use_count(0)
{
  int rval = 0;

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

    try { interval_list.push_back(IntPair(start_ind, end_ind)); }
    catch(const bad_alloc &){
      rval = 1; PSEP_GOTO_CLEANUP("Couldn't push back pair, ");
    }
  }

 CLEANUP:
  if(rval){
    cerr << "IntervalSet constructor failed\n";
    interval_list.clear();
  }
}

SetBank::SetBank(vector<int> &best_tour_nodes,
		 vector<int> &_perm) :
  tour_nodes(best_tour_nodes),
  perm(_perm)
{
  HyperGraph::source_setbank.reset(this);
}

const IntervalSet *SetBank::add_or_increment(IntervalSet &newset)
{
  SetHash::iterator find_it = set_bank.find(newset);

  if(find_it != set_bank.end()){
    find_it->add_use();
  } else {
    newset.add_use();
    set_bank.insert(newset);
    find_it = set_bank.find(newset);    
  }

  return &(*find_it);
}

void SetBank::del_or_decrement(IntervalSet &oldset)
{
  SetHash::iterator find_it = set_bank.find(oldset);

  if(find_it != set_bank.end()){
    find_it->del_use();

    if(find_it->use_count == 0)
      set_bank.erase(find_it);
  }
}


