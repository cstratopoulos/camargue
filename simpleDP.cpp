#include "simpleDP.h"
using namespace std;

int PSEP_SimpleDP::in_subtour_poly(bool *result_p){
  int ncount = G_s.node_count, ecount = G_s.edge_count;
  int end0 = 0;
  double cutval = 2;
  *result_p = false;
  
  for(int end1 = 1; end1 < ncount; end1++){
    if(CCcut_mincut_st(ncount, ecount, &support_elist[0], &support_ecap[0],
		       end0, end1, &cutval, (int **) NULL, (int *) NULL)){
      cerr << "Problem in SimpleDP::separate with Concorde st-cut" << endl;
      return 1;
    }

    if(cutval < 2)
      return 0;
  }

  *result_p = true;
  return 0;
}

void PSEP_SimpleDP::test_build_collection(){
  candidates.build_collection();
  int ncount = G_s.node_count;

  cout << "Printing light teeth...\n";
  for(int i = 0; i < ncount; i++){
    cout << "|||||| i = " << i << ", root = best_tour_nodes[" << i << "] = "
	 << candidates.best_tour_nodes[i] << " |||||\n";
    if(!candidates.light_teeth[i].empty())
      for(list<std::unique_ptr<PSEP_CandTooth::SimpleTooth> >::iterator
	    it = candidates.light_teeth[i].begin();
	  it != candidates.light_teeth[i].end(); it++)
	(*it)->print();
  }
}
