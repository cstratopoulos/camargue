/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief FAST BLOSSOM SEPARATION HEURISTICS
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PSEP_FASTBLOSSOM_H
#define PSEP_FASTBLOSSOM_H

#include "cuts.hpp"
#include "Graph.hpp"

namespace PSEP {

/** Class for primal separation of blossoms via fast standard heuristics.
 * This class performs heuristic primal separation of blossom cuts essentially
 * by using the standard separation heuristics of Padberg-Hong and 
 * Grotschel-Holland.
 */
template<> class Cut<fastblossom> {
public:
  Cut<fastblossom>(std::vector<int> &_delta, std::vector<int> &_edge_marks,
		   Graph &_graph,
		   std::vector<int> &_best_tour_edges,
		   std::vector<double> &_m_lp_edges,
		   std::vector<int> &_support_indices,
		   std::vector<int> &_support_elist,
		   std::vector<double> &_support_ecap,
		   CutQueue<HyperGraph> &fastblossom_queue) :
    delta(_delta), edge_marks(_edge_marks), m_graph(_graph),
    best_tour_edges(_best_tour_edges),
    m_lp_edges(_m_lp_edges),
    support_indices(_support_indices),
    support_elist(_support_elist), support_ecap(_support_ecap),
    local_q(fastblossom_queue.q_capacity),
    external_q(fastblossom_queue) {}

  int cutcall();

protected:
  int separate();
  int add_cuts();

private:
  int build_hypergraph(const fastblossom &blossom_cut);

  int oc_sep(); /**< The Padberg-Hong odd component heuristic. */
  int GH_sep(); /**< The Grotschel-Holland epsilon component heuristic. */

  std::vector<int> &delta;
  std::vector<int> &edge_marks;
  Graph &m_graph;
  
  std::vector<int> &best_tour_edges;
  
  std::vector<double> &m_lp_edges;
  
  std::vector<int> &support_indices;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  CutQueue<fastblossom> local_q;
  CutQueue<HyperGraph> &external_q;
};


}

#endif
