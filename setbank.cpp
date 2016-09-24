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
  
void hashing_test(){
  IntervalSet set1({IntPair(1, 2), IntPair(3,4)});
  IntervalSet set2({IntPair(3,5), IntPair(7,9)});
  IntervalSet dup1({IntPair(3,4), IntPair(1,2)});
  dup1.use_count++;

  std::cout << "Created set 1:\n"; set1.print();
  std::cout << "set1 use count: " << set1.use_count << ", dup1: "
	    << dup1.use_count << "\n";
  std::cout << "Test set1 == dup1: " << (set1 == dup1) << "\n";

  SetHash set_bank;
  set_bank.insert(set1);
  set_bank.insert(set2);
  std::cout << "After adding set1 and 2, set_bank size: "
	    << set_bank.size() << "\n";
  std::cout << "Trying to find dup1: " << set_bank.count(dup1) << "\n";

  sort(dup1.interval_list.begin(), dup1.interval_list.end());

  std::cout << "After sorting, dup1 == set1: " << (set1 == dup1) << "\n";
  std::cout << "Trying to find dup1: " << set_bank.count(dup1) << "\n";
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

    set_bank.insert(newset);
    cout << "Size of set bank is now: " << set_bank.size() << "\n";

    cout << "Continue? y for yes: ";
    cin >> answer;
  }
  cout << "Done testing.\n";  
}

}

using namespace PSEP;

IntervalSet::IntervalSet(vector<int> &nodelist,
			 const vector<int> &best_tour_nodes,
			 const vector<int> &perm) : use_count(0) {
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
