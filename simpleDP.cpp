#include "simpleDP.h"
using namespace std;

void PSEP_SimpleDP::test_build_collection(){
  candidates.build_collection();
  int ncount = G_s.node_count;

  
  cout << "Printing light teeth...\n";
  for(int i = 0; i < ncount; i++){
    cout << "|||||| i = " << i << ", root = best_tour_nodes[" << i << "] = "
	 << candidates.best_tour_nodes[i] << " |||||\n";
    cout << candidates.light_teeth[i].size() << " teeth with this root\n";
    if(!candidates.light_teeth[i].empty())
      for(list<std::unique_ptr<PSEP_CandTooth::SimpleTooth> >::iterator
	    it = candidates.light_teeth[i].begin();
	  it != candidates.light_teeth[i].end(); it++)
	(*it)->print();
  }
}
