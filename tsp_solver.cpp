#include "tsp_solver.h"

using namespace std;

TSP_Solver::TSP_Solver(Graph &graph, const vector<int> &lk_node_indices,
		       PSEP_LP_Prefs _prefs) :
  m_graph(graph),
  best_tour_nodes(lk_node_indices),
  cutcall(m_graph.edges, delta, edge_marks, G_s, best_tour_nodes, m_lp,
	  best_tour_edges, m_lp_edges, support_indices, support_elist,
	  support_ecap, perm, m_graph.edge_lookup),
  LPcore(m_lp, m_graph, m_lp_edges, G_s, support_indices, support_elist,
	 support_ecap, best_tour_edges,best_tour_nodes, perm, m_min_tour_value,
	 island, delta, edge_marks, _prefs),
  print(best_tour_nodes, best_tour_edges, m_lp_edges, m_graph.edges)
{
    //Build the basic LP
    PSEPlp_init (&m_lp);
    PSEPlp_create (&m_lp, "subtour");

    
    if(_prefs.switching_choice == LP::PRICING::SWITCHING::START){
      cout << "Immediate: ";
      LPcore.change_pricing();
    }

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
	LPcore.old_colstat[i] = CPX_BASIC;
	if(m_graph.node_count %2 == 1)
	  continue;	
      }

      if(m_graph.node_count % 2 == 0){
	if(ind0 == m_graph.node_count - 2 && ind1 == m_graph.node_count - 1){
	  LPcore.old_colstat[i] = CPX_AT_UPPER;
	  continue;
	}

	if(ind0 == 0 && ind1 == m_graph.node_count - 2){
	  LPcore.old_colstat[i] = CPX_BASIC;
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

int TSP_Solver::pure_cut(){
  int rval = 0;

  int stat;
  int num_seg, num_2match, num_dp, num_bad, total_cuts = 0;
  int segval, matchval, dpval;
  double segtime, matchtime, dptime;
  int rounds = 0, augrounds = 0;
  bool in_subtour;

  double total_segtime = 0, total_2mtime = 0;
  int total_segcalls = 0, total_2mcalls = 0;

  rval = LPcore.basis_init();
  if(rval) goto CLEANUP;

  cout << "Pivoting until optimality or no more cuts" << endl;

  while(true){
    rounds++;
    augrounds++;
    in_subtour = false;
    
    rval = LPcore.pivot_until_change(&stat);
    if(rval) goto CLEANUP;

    print.pivot(stat);

    if(stat == PIVOT::FATHOMED_TOUR)
      break;

    if(stat == PIVOT::TOUR){
      augrounds = 0;
      rval = LPcore.update_best_tour();
      if(rval)
	goto CLEANUP;
      else {
	cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "!!!AUGMENTED TOUR!!!!" << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "~Call to delete slack cuts should go here~" << endl;
	continue;
      }
    }

    rval = LPcore.pivot_back();
    if(rval) goto CLEANUP;

    segtime = PSEP_zeit();
    segval = cutcall.segment(&num_seg);
    if(segval == 1)
      break;
    segtime = PSEP_zeit() - segtime;
    total_segtime += segtime; total_segcalls++;

    matchtime = PSEP_zeit();
    matchval = cutcall.blossom(2 - num_seg, &num_2match);
    if(matchval == 1)
      break;
    matchtime = PSEP_zeit() - matchtime;
    total_2mtime += matchtime; total_2mcalls++;

    cout << "Added " << num_seg << " segments"
	 << " (in " << segtime << "s)"
	 <<", " << num_2match
	 << " blossoms"
	 << " (in " << matchtime << "s)" << endl;

    dpval = 2; num_dp = 0;
    if(augrounds >= LPcore.prefs.dp_threshold){
      if(num_seg == 0 && stat != PIVOT::SUBTOUR){
	rval = cutcall.in_subtour_poly(&in_subtour);
	if(rval) goto CLEANUP;

	if(in_subtour){
	  dptime = PSEP_zeit();
	  dpval = cutcall.simpleDP(2 - num_2match, &num_dp, &num_bad);
	  if(dpval == 1)
	    break;
	  dptime = PSEP_zeit() - dptime;
	  cout << "*** Added " << num_dp << " simple DP cuts "
	       << "(in " << dptime << "s) (" << num_bad << " bad cuts found)"
	       << " ***\n";
	}
      }
    }

    total_cuts += num_seg + num_2match + num_dp;

    if(num_seg == 0 && num_2match == 0 && num_dp == 0)
      break;
  }

  if(stat != 3)
    cout << "Terminated due to lack of cutting planes after "
	 << rounds << " rounds of separation" << endl;
  cout << total_cuts << " cutting planes added over "
       << rounds << " rounds of separation" << endl;

  cout << "Average time per segcall: "
       << ((double) (total_segtime / total_segcalls)) << "\n";
  cout << "Average time per blossom call: "
       << ((double) (total_2mtime / total_2mcalls)) << "\n";

 CLEANUP:
  if(rval)
    cerr << "Error entry point: TSP_Solver::pure_cut()\n";
  return rval;
}

int TSP_Solver::simple_test(){
  if(LPcore.basis_init())
    return 1;

  int stat, x = 250, num_dp = 0, num_bad = 0;
  int num_seg = 0, segval = 0;
  int num_2match = 0, matchval = 0;
  double segtime, matchtime, routine_start;
  int rounds = 0;
  int total_cuts = 0;

  bool in_sep = false;

  cout << "Pivoting until solution in subtour polytope...\n";

  while(!in_sep){
    if(LPcore.pivot_until_change(&stat))
     return 1;

    print.pivot(stat);

    if(stat == PIVOT::FATHOMED_TOUR){
      cout << "Somehow solved problem, exiting\n";
      return 0;
    }

    if(stat == PIVOT::TOUR){
      if(LPcore.update_best_tour())
	return 1;
      else
	continue;
    }

    if(LPcore.pivot_back())
      return 1;

    segtime = PSEP_zeit();
    segval = cutcall.segment(&num_seg);
    if(segval == 1)
      return 1;
    segtime = PSEP_zeit() - segtime;

    matchtime = PSEP_zeit();
    matchval = cutcall.blossom(250 - num_seg, &num_2match);
    if(matchval == 1)
      return 1;
    matchtime = PSEP_zeit() - matchtime;

   
    cout << "Added " << num_seg << " segments"
	 << " (in " << segtime << "s)"
	 <<", " << num_2match
	 << " blossoms"
	 << " (in " << matchtime << "s)" << endl;

    total_cuts += num_seg + num_2match;

    if(num_seg == 0 && stat != PIVOT::SUBTOUR){
      if(cutcall.in_subtour_poly(&in_sep))
	return 1;
    }

    if(segval + matchval == 4 && !in_sep){
      cout << "No more cuts to add and still not in subtour polytope :'(\n";
      return 1;
    }
    rounds++;
  }

  cout << "Added " << total_cuts << " cuts in "
       << rounds << " rounds\n";

  // print.best_tour_nodes();
  // print.lp_edges();

  if(in_sep){
    cout << "Solution is in subtour polytope, building collection...\n";
    routine_start = PSEP_zeit();
    cutcall.simpleDP(x, &num_dp, &num_bad);
    cout << num_dp << " dp inequalities added\n";
    cout << num_bad << " bad inequalities found\n";
    cout << (PSEP_zeit() - routine_start) << "s finding candidate teeth "
	 << "and building light cutgraph/GH tree\n";
  }

  return 0;
  
}
