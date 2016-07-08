#include "tsp_solver.h"

using namespace std;

TSP_Solver::TSP_Solver(Graph &graph, const vector<int> &lk_node_indices) :
  m_graph(graph),
  best_tour_nodes(lk_node_indices),
  segments(G_s, edge_marks, best_tour_nodes, m_lp),
  blossoms(best_tour_edges, support_indices, support_elist, support_ecap, m_lp)
{
    //Build the basic LP
    PSEPlp_init (&m_lp);
    PSEPlp_create (&m_lp, "subtour");

    /* Build a row for each degree equation */
    for(int i = 0; i < graph.node_count; i++) {
        PSEPlp_new_row (&m_lp, 'E', 2.0);
    }

    /* Build a column for each edge of the Graph */
    int cmatbeg = 0, num_vars = 1, num_non_zero = 2;
    double coefficients[2] = {1.0, 1.0};
    double lower_bound = 0.0;
    double upper_bound = 1.0;
    for(int j = 0; j < graph.edge_count; j++) {
        int *nodes = (int*)graph.edges[j].end;
        double objective_val = (double)graph.edges[j].len;
        PSEPlp_addcols (&m_lp, num_vars, num_non_zero, &objective_val,
                         &cmatbeg, nodes, coefficients, &lower_bound,
			&upper_bound);
    }

    bases.resize(2);
    bases[0].resize(m_graph.edge_count, CPX_AT_LOWER);
    bases[1].resize(m_graph.edge_count, CPX_AT_LOWER);
    old_base = &bases[0][0];
    new_base = &bases[1][0];

    basis_headers.resize(2);
    basis_headers[0].resize(m_graph.node_count);
    basis_headers[1].resize(m_graph.node_count);
    old_header = &basis_headers[0][0];
    new_header = &basis_headers[1][0];
    

    perm.resize(m_graph.node_count);
    for(int i = 0; i < m_graph.node_count; i++){
      perm[best_tour_nodes[i]] = i;
    }

    best_tour_edges.resize(m_graph.edge_count, 0);


    m_min_tour_value = 0;
    for(int i = 0; i < m_graph.edge_count; i++){
      Edge e = m_graph.edges[i];
      int ind0, ind1;
      if(perm[e.end[0]] < perm[e.end[1]]){
	ind0 = perm[e.end[0]]; ind1 = perm[e.end[1]];
      } else {
	ind1 = perm[e.end[0]]; ind0 = perm[e.end[1]];
      }
      
      if(ind1 - ind0 == 1 || (ind0 == 0 && ind1 == m_graph.node_count - 1)){
	best_tour_edges[i] = 1;
	m_min_tour_value += e.len;
	old_base[i] = CPX_BASIC;
	if(m_graph.node_count %2 == 1)
	  continue;	
      }

      if(m_graph.node_count % 2 == 0){
	if(ind0 == m_graph.node_count - 2 && ind1 == m_graph.node_count - 1){
	  old_base[i] = CPX_AT_UPPER;
	  continue;
	}

	if(ind0 == 0 && ind1 == m_graph.node_count - 2){
	  old_base[i] = CPX_BASIC;
	}
      } 
    }

    //Moving some stuff for the DFS up here to save malloc calls
    island.resize(m_graph.node_count);
    delta.resize(m_graph.edge_count, 0);
    edge_marks.resize(m_graph.node_count, 0);
}

TSP_Solver::~TSP_Solver() {
    PSEPlp_free (&m_lp);
}
