#include "n_tooth.hpp"

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <iostream>
#include <algorithm>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::pair;

namespace PSEP {
namespace nu {

CandidateTeeth::CandidateTeeth(PSEP::Data::GraphGroup &_graph_dat,
			       PSEP::Data::BestGroup &_best_dat,
			       PSEP::Data::SupportGroup &_supp_dat) :
  light_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  left_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  right_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  dist_teeth(std::vector<vector<SimpleTooth::Ptr>>(_supp_dat.G_s.node_count)),
  adj_zones(vector<vector<int>>(_supp_dat.G_s.node_count,
				vector<int>(_supp_dat.G_s.node_count, 0))),
  endmark(vector<int>(_supp_dat.G_s.node_count, CC_LINSUB_BOTH_END)),
  graph_dat(_graph_dat),
  best_dat(_best_dat),
  supp_dat(_supp_dat)
{
  PSEP::SupportGraph &G_s = supp_dat.G_s;
  int ncount = G_s.node_count;
  vector<int> &perm = best_dat.perm;
  vector<int> &tour = best_dat.best_tour_nodes;

  for(int root_ind = 0; root_ind < ncount; ++root_ind){
    int actual_vx = tour[root_ind];
    PSEP::SNode x = G_s.nodelist[actual_vx];
    
    for(int k = 0; k < x.s_degree; ++k){
      int end_ind = perm[x.adj_objs[k].other_end];
      
      adj_zones[root_ind][end_ind] = 1;
    }

    int label = 0;
    for(int &i : adj_zones[root_ind]){
      if(i == 1) ++label;
      i = label;
    }
  }
}
  
  

}
}
