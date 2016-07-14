#include "simpleDP.h"
using namespace std;

void PSEP_SimpleDP::build_light_cuttree(){
  int num_cutnodes = 0;
  int current_index, star_index;
  bool found_parent = false;
  list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::iterator parent;
  list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::iterator child;

  cut_elist.clear();
  cut_ecap.clear();
  light_nodes.clear();
  cut_marks.clear();
  node_marks.clear();

  //NULL shared ptrs for the degree eqn on each root
  light_nodes.resize(candidates.light_teeth.size());
  cout << "Use count of one of the null degree eqn nodes: "
       << light_nodes.begin()->use_count() << " = first one, "
       << light_nodes.end()->use_count() << " = last one\n";

  current_index = light_nodes.size();

  cout << "Checking the use count on all the shared ptrs...";
  for(int i = 0; i < candidates.light_teeth.size(); i++)
    for(list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::iterator it
	  = candidates.light_teeth[i].begin();
	it != candidates.light_teeth[i].end(); it++)
      if(it->use_count() != 1)
	cout << "\nFound one w wrong use count???\n";
  cout << "Done\n";

  //number a light node for each tooth ineq, add to vector
  cout << "Building non-NULL list of light nodes, checking use counts....";
  for(int i = 0; i < candidates.light_teeth.size(); i++){
    for(list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::iterator
	  it = candidates.light_teeth[i].begin();
	it != candidates.light_teeth[i].end(); it++){
      light_nodes.push_back(*it);
      (*it)->node_index = current_index++;
      if(it->use_count() != 2)
	cout << "\n found one with wrong use count????\n";
    }
  }
  cout << "Done\n";

  light_nodes.emplace_back((PSEP_CandTooth::SimpleTooth *) NULL);
  num_cutnodes = light_nodes.size();
  star_index = num_cutnodes - 1; // the special central node
  
  node_marks.resize(num_cutnodes, false);

  //add edges between every root and the central node, the degree eqns
  for(int i = 0; i < candidates.light_teeth.size(); i++){
    cut_elist.push_back(i);
    cut_elist.push_back(star_index);
    cut_ecap.push_back(0);
  }

  //add simple tooth inequality edges
  for(int i = 0; i < candidates.light_teeth.size(); i++){
    if(candidates.light_teeth[i].empty()) continue;

    child = candidates.light_teeth[i].begin();
    cut_elist.push_back((*child)->node_index);
    cut_elist.push_back(i);
    cut_ecap.push_back((*child)->slack);
    node_marks[i] = !(node_marks[i]);
    node_marks[(*child)->node_index] = !(node_marks[(*child)->node_index]);
    //bits in node_marks are complemented to count odd or evenness

    child++; //for loop initialized one after begin
    for(; child != candidates.light_teeth[i].end(); child++){

      parent = child; //while loop initializer
      while(true){
	parent --;

	found_parent =
	  PSEP_CandTooth::SimpleTooth::C_body_subset(**child, **parent);

	if(found_parent){//edge between smallest parent
	  cut_elist.push_back((*child)->node_index);
	  cut_elist.push_back((*parent)->node_index);
	  node_marks[(*child)->node_index]
	    = !(node_marks[(*child)->node_index]);
	  node_marks[(*parent)->node_index]
	    = !(node_marks[(*parent)->node_index]);
	  cut_ecap.push_back((*child)->slack);
	  break;
	}

	if(parent == candidates.light_teeth[i].begin()) break;
      }

      if(!found_parent){//edge is between degree node if no parent found
	cut_elist.push_back((*child)->node_index);
	cut_elist.push_back(i);
	node_marks[i] = !(node_marks[i]);
	node_marks[(*child)->node_index] = !(node_marks[(*child)->node_index]);
	cut_ecap.push_back((*child)->slack);
      }
    }
  }

  //setting vector of marks for use by concorde gomoryhu
  for(int i = 0; i < num_cutnodes; i++)
    if(node_marks[i])
      cut_marks.push_back(i);
  
}

void PSEP_SimpleDP::test_build_collection(){
  candidates.build_collection();
  int ncount = G_s.node_count;
  int light_total = 0, heavy_total = 0;

  cout << "Printing light teeth...\n";
  for(int i = 0; i < ncount; i++){
    light_total += candidates.light_teeth[i].size();
  }

  cout << light_total << " light teeth in total\n";

  cout << "Printing heavy teeth...\n";
  for(int i = 0; i < ncount; i++){
    heavy_total += candidates.heavy_teeth[i].size();
  }

  cout << heavy_total << " heavy teeth in total\n";

  cout << "Total number of teeth with slack less than one: "
       << light_total + heavy_total << endl;

  cout << "Ratio of light teeth to n^2 is : "
       << ((double) light_total / (ncount * ncount)) << "\n";

  cout << "Ratio of heavy teeth to n^3 is : "
       << ((double) heavy_total / (ncount * ncount * ncount)) << "\n";

  cout << "Ratio of total number of teeth to n^3 is: "
       << ((double) (light_total + heavy_total) / (ncount * ncount * ncount))
       << "\n";

  cout << "Now calling build light cuttree..................\n";
  build_light_cuttree();

}
