#ifndef PSEP_LPWORK_H
#define PSEP_LPWORK_H

#include "lp.h"
#include "Graph.h"

class PSEP_LP_Core {
 public:
 PSEP_LP_Core(PSEPlp & _lp, Graph & _graph, std::vector<double> & _lp_edges,
	      SupportGraph & _G, std::vector<int> & _sup_inds,
	      std::vector<int> & _sup_elist, std::vector<double> & _sup_ecap,
	      std::vector<int> & _tour_edges, std::vector<int> & _tour_nodes,
	      std::vector<int> & _perm,
	      double & _tour_val, std::vector<int> & _island,
	      std::vector<int> & _delta, std::vector<int> & _edge_marks,
	      PSEP_LP_Prefs _prefs) :
  m_lp(_lp), m_graph(_graph), prefs(_prefs), m_lp_edges(_lp_edges), G_s(_G),
    support_indices(_sup_inds), support_elist(_sup_elist),
    support_ecap(_sup_ecap), best_tour_edges(_tour_edges),
    best_tour_nodes(_tour_nodes), perm(_perm), m_min_tour_value(_tour_val),
    island(_island), delta(_delta), edge_marks(_edge_marks)
    {
    old_colstat.resize(m_graph.edge_count, CPX_AT_LOWER);
    old_rowstat.resize(m_graph.node_count, CPX_AT_LOWER);      
    }
    
  bool is_dual_feas();
  bool is_integral();

  int pivot();
  int pivot_back();

  int basis_init();
  
  int pivot_until_change(int *pivot_status_p);


  double get_obj_val();
  double set_edges();
  double set_support_graph();

  int update_best_tour();

  int prune_cuts(int *num_removed);

  void change_pricing();

 private:
  PSEPlp &m_lp;
  Graph &m_graph;

  PSEP_LP_Prefs prefs;

  friend class TSP_Solver;
  friend class PSEP_Cutcall;
  std::vector<int> old_colstat;
  std::vector<int> old_rowstat;

  std::vector<double> &m_lp_edges;

  SupportGraph &G_s;
  
  std::vector<int> &support_indices;
  std::vector<int> &support_elist;
  std::vector<double> &support_ecap;

  std::vector<int> &best_tour_edges;
  std::vector<int> &best_tour_nodes;
  std::vector<int> &perm;

  double &m_min_tour_value;

  std::vector<int> &island;
  std::vector<int> &delta;
  std::vector<int> &edge_marks;
  
};


#endif
