#ifndef CMR_WITNESS_H
#define CMR_WITNESS_H

#include "datagroups.hpp"
#include "tooth.hpp"
#include "process_cuts.hpp"


#include <vector>

extern "C" {
#include <concorde/INCLUDE/cut.h>
}

namespace CMR {
namespace Sep {

/** Class for building miniature simple DP witness cutgraphs. */
class DPwitness {
public:
  /** Construct a mini cutgraph induced by a partition. */
  DPwitness(CandidateTeeth &cands,
            const std::vector<int> &partition_nodes);
  ~DPwitness();
  
  /** Create a cutgraph and grab odd cuts from it. */
  bool simple_DP_sep(CutQueue<dominoparity> &domino_q);
  
private:
  void build_light_tree();  /**< Build the tooth inequality tree. */
  void add_web_edges(); /**< Add nonnegativity inequality edges. */
  void build_gh_tree(); /**< Construct the Gomory-Hu tree. */

  void dfs_odd_cuts(CC_GHnode *n); /**< Search tree for odd cuts. */

  /** Expand a Gomory-Hu node into a cut. */
  void expand_cut(CC_GHnode *n, std::vector<int> &cut_nodes);

  /** Get simple DP inequalities from fundamental cuts of the GH tree. */
  void grab_dominos(CutQueue<dominoparity> &domino_q);

  std::vector<std::vector<SimpleTooth>> light_teeth;

  const Data::SupportGroup &supp_dat;

  const std::vector<int> &perm;

  std::vector<SimpleTooth*> cutgraph_nodes;

  std::vector<int> cut_elist;
  std::vector<double> cut_ecap;
  std::vector<int> odd_nodes_list;
  std::vector<bool> node_marks;

  std::vector<int> cutgraph_delta;
  std::vector<int> cg_delta_marks;

  CutQueue<CC_GHnode *> CC_gh_q;

  CC_GHtree gh_tree;
  
  void grab_cuts(CutQueue<dominoparity> &domino_q);

};

}
}

#endif
