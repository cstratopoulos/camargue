/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief BUILDING SIMPLE DP WITNESS CUTGRAPH
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_DPGRAPH_HPP
#define PSEP_DPGRAPH_HPP

#include "tooth.hpp"
#include "cuts.hpp"

#include <vector>
#include <string>
#include <fstream>

extern "C" {
#include <concorde/INCLUDE/cut.h>
}


/** @def PSEP_DO_VIZ
 * @brief Macro for conditional compilation of graphviz DOT files.
 * If defined, the DPCutGraph constructor must provide an output file name.
 * This will be used to open an output file stream, which will contain a 
 * graphviz (.gv) file describing the light cut tree and full witness cutgraph.
 * This can then be processed with (suggested) fdp, circo, etc. 
 * @remark Comments in the gv file will indicate the beginning of the 
 * nonnegativity edges. You may wish to delete these entirely; with them 
 * included it is hard to get a graph that is not an unintelligible scribble.
 * @warning This should be undef'd for all but very small examples in Catch
 * test cases. 
 */
#undef PSEP_DO_VIZ

namespace PSEP {

/** Class for building light simple DP witness cutgraphs.
 * As per Fleischer, Letchford, and Lodi (2006), this class will build the 
 * witness graph for detecting light simple DP inequalities. Then, it builds
 * a Gomory-Hu cut tree from the witness graph, from which an odd cut of weight
 * less than one corresponds to a violated light simple DP inequality.
 */
class DPCutGraph {
public:
  DPCutGraph(
#ifdef PSEP_DO_VIZ
	     std::string _ofname,
#endif
	     std::vector<std::vector<PSEP::SimpleTooth::Ptr>> &_teeth,
	     PSEP::CandidateTeeth &_cands,
	     const std::vector<int> &_perm,
	     const SupportGraph &_G_s,
	     const std::vector<int> &_support_elist,
	     const std::vector<double> &_support_ecap);
  ~DPCutGraph();

#ifdef PSEP_DO_VIZ
  std::string ofname;
  std::ofstream cg_out;
#endif

  int simple_DP_sep(PSEP::CutQueue<PSEP::dominoparity> &domino_q);

  //  private:
  int grab_cuts(PSEP::CutQueue<PSEP::dominoparity> &domino_q);

  int build_light_tree();
  int add_web_edges();
  int call_concorde_gomoryhu();

  void dfs_odd_cuts(CC_GHnode *n);
  void expand_cut(CC_GHnode *n, std::vector<int> &cut_shore_nodes);
  
  const std::vector<std::vector<PSEP::SimpleTooth::Ptr>> &light_teeth;
  PSEP::CandidateTeeth &cands;

  const SupportGraph &G_s;
  const std::vector<int> &support_elist;
  const std::vector<double> &support_ecap;
  
  const std::vector<int> &perm;
  
  std::vector<PSEP::SimpleTooth*> cutgraph_nodes;
  
  std::vector<int> cut_elist;
  std::vector<double> cut_ecap;
  std::vector<int> odd_nodes_list;
  std::vector<bool> node_marks;

  std::vector<int> cutgraph_delta;
  std::vector<int> cg_delta_marks;

  PSEP::CutQueue<CC_GHnode *> CC_gh_q;

  CC_GHtree gh_tree;

  
};

}

#endif
