/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief EXACT PRIMAL SUBTOUR SEPARATION
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_SEGCUTS_H
#define CMR_SEGCUTS_H

#include "cuts.hpp"
#include "Graph.hpp"

#include<vector>

namespace CMR {

/** Class for exact primal subtour separation.
 * Can be instantiated with data about a best tour and an lp solution, and used
 * to perform primal subtour separation using the segments algorithm of
 * Applegate et al.
 */
template<> class Cut<seg> {
public:
  Cut<seg>(std::vector<int> &_best_tour_nodes, std::vector<int> &_perm,
	   SupportGraph &_G_s, std::vector<int> &_support_elist,
	   std::vector<double> &_support_ecap,
	   CutQueue<HyperGraph> &segment_queue):
    best_tour_nodes(_best_tour_nodes), perm(_perm),
    G_s(_G_s), support_elist(_support_elist), support_ecap(_support_ecap),
    local_q(segment_queue.q_capacity),
    external_q(segment_queue){}

  int cutcall();

protected:
  int separate();
  int add_cuts();
    
private:
  int build_hypergraph(const seg& seg_cut);

  /** Callback function for use by Concorde linsub function.
   * This function is the callback to CCcut_linsub in Concorde. For more
   * information see segments.c in Concorde. 
   */
  static int linsub_callback(double cut_val, int cut_start, int cut_end,
			     void *cb_data);

  std::vector<int> &best_tour_nodes;
  std::vector<int> &perm;

  SupportGraph &G_s;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  CutQueue<seg> local_q;
  CutQueue<HyperGraph> &external_q;
};

}

#endif
