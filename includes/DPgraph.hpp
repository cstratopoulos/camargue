#ifndef PSEP_DPGRAPH_HPP
#define PSEP_DPGRAPH_HPP

#include "tooth.hpp"
#include "cuts.hpp"

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

#include <vector>

namespace PSEP {

class DPCutGraph {
public:
  DPCutGraph(std::vector<std::vector<PSEP::SimpleTooth::Ptr>> &_teeth,
	     const std::vector<int> &_perm,
	     const SupportGraph &_G_s);
  ~DPCutGraph();

  int grab_cuts(PSEP::CutQueue<PSEP::dominoparity> &domino_q);

  //private:
  int build_light_tree();
  int add_web_edges();
  int call_concorde_gomoryhu();

  void dfs_odd_cuts(CC_GHnode *n, int &cutcount);
  
  const std::vector<std::vector<PSEP::SimpleTooth::Ptr>> &light_teeth;

  const SupportGraph &G_s;
  const std::vector<int> &perm;
  
  std::vector<PSEP::SimpleTooth*> cutgraph_nodes;
  
  std::vector<int> cut_elist;
  std::vector<double> cut_ecap;
  std::vector<int> odd_nodes_list;
  std::vector<bool> node_marks;

  CC_GHtree gh_tree;

  
};

}

#endif
