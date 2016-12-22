/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief BUILDING SIMPLE DP WITNESS CUTGRAPH
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_DPGRAPH_HPP
#define CMR_DPGRAPH_HPP

#include "tooth.hpp"
#include "process_cuts.hpp"

#include <vector>
#include <string>
#include <fstream>

extern "C" {
#include <concorde/INCLUDE/cut.h>
}


/** @def CMR_DO_VIZ
 * @brief Macro for conditional compilation of graphviz DOT files.
 * If defined, the DPCutGraph constructor must provide an output file name.
 * This will be used to open an output file stream, which will contain a 
 * graphviz (.gv) file describing the light cut tree.
 * This can then be processed with (suggested) fdp, circo, etc. 
 * @warning This should be undef'd for all but very small examples in Catch
 * test cases. 
 */
#undef CMR_DO_VIZ

namespace CMR {

/** Class for building light simple DP witness cutgraphs.
 * As per Fleischer, Letchford, and Lodi (2006), this class will build the 
 * witness graph for detecting light simple DP inequalities. Then, it builds
 * a Gomory-Hu cut tree from the witness graph, from which an odd cut of weight
 * less than one corresponds to a violated light simple DP inequality.
 */
class DPCutGraph {
public:
  DPCutGraph(
#ifdef CMR_DO_VIZ
	     std::string _ofname,
#endif
	     CMR::CandidateTeeth &_cands);
  ~DPCutGraph();

#ifdef CMR_DO_VIZ
  std::string ofname;
  std::ofstream cg_out;
#endif

  int simple_DP_sep(Sep::CutQueue<Sep::dominoparity> &domino_q);

  //  private:
  int grab_cuts(Sep::CutQueue<Sep::dominoparity> &domino_q);

  int build_light_tree();
  int add_web_edges();
  int call_concorde_gomoryhu();

  void dfs_odd_cuts(CC_GHnode *n);
  void expand_cut(CC_GHnode *n, std::vector<int> &cut_shore_nodes);
  
  const std::vector<std::vector<CMR::SimpleTooth::Ptr>> &light_teeth;
  CMR::CandidateTeeth &cands;

  const SupportGraph &G_s;
  const std::vector<int> &support_elist;
  const std::vector<double> &support_ecap;
  
  const std::vector<int> &perm;
  
  std::vector<CMR::SimpleTooth*> cutgraph_nodes;
  
  std::vector<int> cut_elist;
  std::vector<double> cut_ecap;
  std::vector<int> odd_nodes_list;
  std::vector<bool> node_marks;

  std::vector<int> cutgraph_delta;
  std::vector<int> cg_delta_marks;

  Sep::CutQueue<CC_GHnode *> CC_gh_q;

  CC_GHtree gh_tree;

  
};

}

#endif
