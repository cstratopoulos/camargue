#include "simpleDP.h"
using namespace std;

int PSEP_SimpleDP::separate(const int max_cutcount){
  int rval = 0;
  int ncount = G_s.node_count;
  int light_total = 0, heavy_total = 0;

  toothlists.clear();

  candidates.build_collection();

  for(int i = 0; i < ncount; i++){
    light_total += candidates.light_teeth[i].size();
    cout << candidates.light_teeth[i].size() << " light teeth with root "
	 << i << ", ratio / num nodes: "
	 << ((double) candidates.light_teeth[i].size() / ncount) << "\n";
    heavy_total += candidates.heavy_teeth[i].size();
  }

  cout << light_total << " light teeth, " << heavy_total << " heavy teeth\n";
  cout << "Light teeth / n^2 = "
       << ((double) light_total / (ncount * ncount)) << "\n";
  cout << "Total num teeth / n^3 = "
       << ((double) (light_total + heavy_total) / (ncount * ncount * ncount))
       << "\n";

  cout << "Building light cuttree...."; build_light_cuttree();
  cout << "Done\n";

  cout << "Adding web edges...."; add_web_edges(); cout << "Done\n";

  /*
  cout << "Calling concorde and building toothlists...";
  rval = call_CC_gomoryhu(max_cutcount);
  if(rval)
    goto CLEANUP;
  cout << "Done\n";

  cout << "Number of inequalities for consideration: "
       << toothlists.size() << "\n";
  */

 CLEANUP:
  if(rval)
    cerr << "Error entry point: SimpleDP::separate\n";
  return rval;
}

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

  current_index = light_nodes.size();

  //number a light node for each tooth ineq, add to vector
  for(int i = 0; i < candidates.light_teeth.size(); i++){
    for(list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::iterator
	  it = candidates.light_teeth[i].begin();
	it != candidates.light_teeth[i].end(); it++){
      light_nodes.push_back(*it);
      (*it)->node_index = current_index++;
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

void PSEP_SimpleDP::add_web_edges(){
  double lp_weight;
  bool end1_in = false, end0_in = false;

  for(int end0 = 0; end0 < G_s.node_count; end0++){
    for(int j = 0; j < G_s.nodelist[end0].s_degree; j++){
      end0_in = false; end1_in = false;
      int end1 = G_s.nodelist[end0].adj_objs[j].other_end;
      if(end0 > end1) continue;
      lp_weight = G_s.nodelist[end0].adj_objs[j].lp_weight;
      if(!candidates.light_teeth[end0].empty()){
	for(list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::reverse_iterator
	      root_end0 = candidates.light_teeth[end0].rbegin();
	    root_end0 != candidates.light_teeth[end0].rend(); root_end0++){
	  end1_in = (*root_end0)->body_contains(perm[end1]);
	  if(end1_in){
	    cut_elist.push_back((*root_end0)->node_index);
	    break;
	  }
	}
      }
      if(!end1_in)
	cut_elist.push_back(end0);

      if(!candidates.light_teeth[end1].empty()){
	for(list<shared_ptr<PSEP_CandTooth::SimpleTooth> >::reverse_iterator
	      root_end1 = candidates.light_teeth[end1].rbegin();
	    root_end1 != candidates.light_teeth[end1].rend(); root_end1++){
	  end0_in = (*root_end1)->body_contains(perm[end0]);
	  if(end0_in){
	    cut_elist.push_back((*root_end1)->node_index);
	    break;
	  }
	}
      }
      if(!end0_in)
	cut_elist.push_back(end1);

      cut_ecap.push_back(lp_weight);
    }
  }
}

int PSEP_SimpleDP::call_CC_gomoryhu(const int max_cutcount){
  int rval = 0;
  int ncount = light_nodes.size();
  int ecount = cut_ecap.size();
  int markcount = cut_marks.size();
  
  CCrandstate rstate;
  CC_GHtree T;
  int seed = (int) PSEP_real_zeit();
  CCutil_sprand(seed, &rstate);
  CCcut_GHtreeinit(&T);

  CC::GH::cut_pq node_pq;

  toothlists.clear();

  if(cut_elist.empty() || cut_ecap.empty() || cut_marks.empty()){
    cerr << "Passed empty vector to gomoryhu\n";
    rval = 1; goto CLEANUP;
  }

  rval = CCcut_gomory_hu(&T, ncount, ecount, &cut_elist[0], &cut_ecap[0],
			 markcount, &cut_marks[0], &rstate);

  if(rval){
    cerr << "Problem calling CCcut_gomory_hu\n";
    goto CLEANUP;
  }

  CC::GH::get_all_toothlists(&T, max_cutcount, node_pq, toothlists);

 CLEANUP:
  CCcut_GHtreefree(&T);
  if(rval)
    cerr << "Error entry point: SimpleDP::call_CC_gomoryhu\n";
  return rval;
}
