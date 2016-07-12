#ifndef PSEP_LPWORK_H
#define PSEP_LPWORK_H

#include "lp.h"
#include "Graph.h"

class PSEP_LP_Core {
 public:
 PSEP_LP_Core(PSEPlp & _lp, Graph & _graph, std::vector<int> & _sup_inds,
	      std::vector<int> & _sup_elist, std::vector<int> _sup_ecap,
	      std::vector<int> & _tour_edges, std::vector<int> & _tour_nodes,
	      std::vector<int> & _perm,
	      double & _tour_val, std::vector<int> & _island,
	      std::vector<int> & _delta, std::vector<int> & _edge_marks) :
    m_lp(_lp), m_graph(_graph),
    support_indices(_sup_inds), support_elist(_sup_elist),
    support_ecap(_sup_ecap), best_tour_edges(_tour_edges),
    best_tour_nodes(_tour_nodes), perm(_perm), m_min_tour_val(_tour_val)
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

  void enable_devex();

 private:
  PSEPlp &m_lp;
  Graph &m_graph;

  std::vector<int> old_colstat;
  std::vector<int> old_rowstat;

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
