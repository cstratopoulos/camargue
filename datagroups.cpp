#include "datagroups.hpp"

#include<algorithm>
#include<unordered_map>
#include<iostream>
#include<iomanip>
#include<vector>

#include<cmath>

extern "C" {
#include <concorde/INCLUDE/linkern.h>
#include <concorde/INCLUDE/util.h>
#include <concorde/INCLUDE/edgegen.h>
}

#include "PSEP_util.hpp"

using namespace std;
using namespace PSEP::Data;

GraphGroup::GraphGroup(const string &fname, RandProb &randprob,
		       unique_ptr<CCdatagroup> &dat,
		       const bool sparse, const int quadnearest){
  int rval = 0;
  CCdatagroup *rawdat = dat.get();
  char *filestring = const_cast<char *>(fname.data());
  int *elist = (int *) NULL;

  CCutil_init_datagroup(rawdat);

  if(!fname.empty()){
    rval = CCutil_gettsplib(filestring, &(m_graph.node_count), rawdat);
    if (rval) PSEP_GOTO_CLEANUP("CCutil_gettsplib failed, ");
  }
  else {
    m_graph.node_count = randprob.nodecount;
    CCrandstate rstate;
    int use_gridsize = randprob.gridsize;
    int allow_dups = 1;
    if(randprob.seed == 0) randprob.seed = (int) real_zeit();
    cout << "Random seed: " << randprob.seed << "\n";
    CCutil_sprand(randprob.seed, &rstate);
    rval = CCutil_getdata((char *) NULL, 1, CC_EUCLIDEAN,
			  &(m_graph.node_count),
			  rawdat, use_gridsize, allow_dups, &rstate);
    if(rval) PSEP_GOTO_CLEANUP("CCutil_getdata randprob failed, ");
  }

  if(!sparse){  
    m_graph.edge_count = (m_graph.node_count * (m_graph.node_count - 1)) / 2;
    try{ m_graph.edges.resize(m_graph.edge_count); }
    catch(const std::bad_alloc &){
      rval = 1; PSEP_GOTO_CLEANUP("Out of memory for m_graph.edges, ");
    }
  
    { int e_index = 0; 
      for (int i = 0; i < m_graph.node_count; i++){
	for (int j = i+1; j < m_graph.node_count; j++){
	  Edge e(i, j, CCutil_dat_edgelen(i, j, rawdat));
	  
	  m_graph.edges[e_index] = e;
	  m_graph.edge_lookup.emplace(IntPair(i,j), e_index);
	  e_index++;
	}
      }}
  } else {
    cout << "    GENERATING SPARSE GRAPH ONLY, ";
    CCedgegengroup plan;
    CCrandstate rstate;
    int edgegen_seed = (randprob.seed == 0) ? ((int) real_zeit()) :
      randprob.seed;
    
    CCutil_sprand(edgegen_seed, &rstate);
    CCedgegen_init_edgegengroup(&plan);
    plan.linkern.count = 10;
    plan.linkern.quadnearest = 5;
    plan.linkern.greedy_start = 0;
    plan.linkern.nkicks = (m_graph.node_count / 100) + 1;
    plan.quadnearest = quadnearest;

    cout << plan.linkern.count << " LK trials, quad "
	 << plan.quadnearest << "-nearest, seed " << edgegen_seed
	 << "\n";

    rval = CCedgegen_edges(&plan, m_graph.node_count, rawdat, NULL,
			   &(m_graph.edge_count), &elist, 1, &rstate);
    if(rval) PSEP_GOTO_CLEANUP("Problem with CCedgegen_edges, ");
    
    try{ m_graph.edges.resize(m_graph.edge_count); }
    catch(const std::bad_alloc &){
      rval = 1; PSEP_GOTO_CLEANUP("Out of memory for m_graph.edges, ");
    }

    for(int i = 0; i < m_graph.edge_count; i++){
      Edge e(elist[2 * i], elist[(2 * i) + 1],
	     CCutil_dat_edgelen(elist[2 * i], elist[(2 * i) + 1], rawdat));
      
      m_graph.edges[i] = e;
      m_graph.edge_lookup.emplace(IntPair(e.end[0], e.end[1]), i);
    }

    cout << "    " << m_graph.edge_count
	 << " edges in sparse graph, ratio to nodes: " << setprecision(2)
	 << ((double) m_graph.edge_count / m_graph.node_count)
	 << setprecision(6) << "\n";
  }

  try{
  island.resize(m_graph.node_count);
  delta.resize(m_graph.edge_count, 0);
  edge_marks.resize(m_graph.node_count, 0);
  } catch(const std::bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for dfs vectors, ");
  }

 CLEANUP:
  if(elist) free(elist);
  if(rval){
    cerr << "Problem in GraphGroup constructor\n";
    m_graph.node_count = 0;
  }
}

BestGroup::BestGroup(Graph &m_graph, vector<int> &delta,
		     unique_ptr<CCdatagroup> &dat, const int user_seed){
  int rval = 0;
  CCrandstate rand_state;
  CCedgegengroup plan;
  CCdatagroup *rawdat = dat.get();
  int ncount = m_graph.node_count;
  int ecount = 0;
  int *elist = (int *) NULL;
  int tcount = 0;
  int *tlist = (int *) NULL;
  int *cyc = (int *) NULL;
  double bestval, val;
  int trials = 4;// (user_seed == 0) ? 4 : 0;
  int silent = 1;
  int kicks = 5 * ncount;
  int istour;
  int seed = (user_seed == 0) ? ((int) real_zeit()) : user_seed;
  bool sparse = (m_graph.edge_count < (ncount * (ncount - 1)) / 2);

  cout << "LK seed: " << seed << ", " << trials << " trials\n";

  bestval = INFINITY;

  //code copies from static int find_tour from concorde
  CCutil_sprand(seed, &rand_state);

  cyc = CC_SAFE_MALLOC(ncount, int); 
  if(!cyc){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for cyc, ");
  }

  try {
  best_tour_nodes.resize(ncount);
  perm.resize(ncount);
  best_tour_edges.resize(m_graph.edge_count, 0);
  } catch(const std::bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for BestGroup vectors, ");
  }

  if(!sparse){
    cout << "Doing non-sparse quad-nearest and greedy tour\n";
    CCedgegen_init_edgegengroup (&plan);
    plan.quadnearest = 2;
    rval = CCedgegen_edges (&plan, ncount, rawdat, (double *) NULL, &ecount,
			    &elist, silent, &rand_state);
    if(rval) PSEP_GOTO_CLEANUP("CCedgegen_edges failed, ");
    plan.quadnearest = 0;

    plan.tour.greedy = 1;
    rval = CCedgegen_edges (&plan, ncount, rawdat, (double *) NULL, &tcount,
			    &tlist, silent, &rand_state);
    if(rval) PSEP_GOTO_CLEANUP("CCedgegen_edges failed, ");

    if (tcount != ncount) {
      rval = 1; PSEP_GOTO_CLEANUP("Wrong edgeset from CCedgegen_edges, ");
    }

    rval = CCutil_edge_to_cycle (ncount, tlist, &istour, cyc);
    if(rval) PSEP_GOTO_CLEANUP("CCutil_edge_to_cycle failed, ");
  
    if (istour == 0) {
      rval = 1; PSEP_GOTO_CLEANUP("Starting tour has an error, ");
    }
    CC_FREE (tlist, int);
  } else {
    ecount = m_graph.edge_count;
    elist = CC_SAFE_MALLOC(2 * ecount, int);
    if(!elist){
      rval = 1; PSEP_GOTO_CLEANUP("Out of memory for elist, ");
    }

    for(int i = 0; i < ecount; i++){
      elist[2 * i] = m_graph.edges[i].end[0];
      elist[(2 * i) + 1] = m_graph.edges[i].end[1];
    }
  }

  //in a sparse graph we do not bother putting a greedy tour as incycle
  rval = CClinkern_tour (ncount, rawdat, ecount, elist, ncount, kicks,
			 sparse ? (int *) NULL : cyc, //
			 &best_tour_nodes[0], &bestval, silent, 0.0, 0.0,
			 (char *) NULL,
			 CC_LK_GEOMETRIC_KICK, &rand_state);
  if(rval) PSEP_GOTO_CLEANUP("CClinkern_tour failed, ");
  
  //end of copied code (from find_tour)
  cout << "LK initial run: " << bestval << endl;

  //begin copied code from find_good_tour in tsp_call.c  
  for(int i = 0; i < trials; i++){
    rval = CClinkern_tour(ncount, rawdat, ecount, elist, ncount, kicks,
			  (int *) NULL, cyc, &val, 1, 0.0, 0.0,
			  (char *) NULL, CC_LK_GEOMETRIC_KICK, &rand_state);
  if(rval) PSEP_GOTO_CLEANUP("CClinkern_tour failed, ");
  
    cout << "LK run " << i << ": " << val << "\n";
    if(val < bestval){
      for(int j = 0; j < ncount; j++)
	best_tour_nodes[j] = cyc[j];
      bestval = val;
    }
  }
  
  if (trials > 0){
    rval = CClinkern_tour(ncount, rawdat, ecount, elist, ncount, 2 * kicks,
			  &best_tour_nodes[0], cyc, &bestval, 1, 0.0, 0.0,
			  (char *) NULL,
			  CC_LK_GEOMETRIC_KICK, &rand_state);
    if(rval) PSEP_GOTO_CLEANUP("CClinkern_tour failed, ");

    cout << "LK run from best tour: " << bestval << "\n";
    for(int j = 0; j < ncount; j++)
      best_tour_nodes[j] = cyc[j];
  }

  for(int i = 0; i < m_graph.node_count; i++)
    perm[best_tour_nodes[i]] = i;
  
  min_tour_value = bestval;

  { int missing = 0;
  for(int i = 0; i < ncount; i++){
    int end0 = fmin(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
    int end1 = fmax(best_tour_nodes[i], best_tour_nodes[(i + 1) % ncount]);
    IntPairMap::const_iterator edge_it =
      m_graph.edge_lookup.find(IntPair(end0, end1));
    
    if(edge_it == m_graph.edge_lookup.end()){
      missing++;

      Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, rawdat));

      m_graph.edges.push_back(e);
      m_graph.edge_lookup[IntPair(end0, end1)] = m_graph.edges.size() - 1;
      m_graph.edge_count += 1;
      best_tour_edges.push_back(0);
      delta.push_back(0);
      edge_it = m_graph.edge_lookup.find(IntPair(end0, end1));
    }

    int edge_index = edge_it->second;    
    best_tour_edges[edge_index] = 1;
  }
  cout << "Added " << missing << " additional edges. ";
  }

  if((ncount % 2) == 0){
    int end0 = fmin(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
    int end1 = fmax(best_tour_nodes[0], best_tour_nodes[ncount - 2]);
    IntPairMap::const_iterator edge_it =
      m_graph.edge_lookup.find(IntPair(end0, end1));

    if(edge_it == m_graph.edge_lookup.end()){
      cout << "Adding extra edge just for basis.";
      Edge e(end0, end1, CCutil_dat_edgelen(end0, end1, rawdat));

      m_graph.edges.push_back(e);
      m_graph.edge_lookup[IntPair(end0, end1)] = m_graph.edges.size() - 1;
      m_graph.edge_count += 1;
      best_tour_edges.push_back(0);
      delta.push_back(0);
    }    
  }
  cout << "\n";

 CLEANUP:
  CC_IFFREE (cyc, int);
  CC_IFFREE (elist, int);
  CC_IFFREE (tlist, int);
  CCutil_freedatagroup(dat.get());
  if(rval){
    cerr << "Problem in BestGroup constructor\n";
    best_tour_edges.clear();
    best_tour_nodes.clear();
    perm.clear();
  }
}

LPGroup::LPGroup(const Graph &m_graph, PSEP::LP::Prefs &_prefs,
			   const vector<int> &perm){
  int rval = 0;
  int cmatbeg = 0, num_vars = 1, num_non_zero = 2;
  double coefficients[2] = {1.0, 1.0};
  double lower_bound = 0.0;
  double upper_bound = 1.0;

  
  //Build the basic LP
  rval = PSEPlp_init (&m_lp);
  if(rval) goto CLEANUP;

  rval = PSEPlp_create (&m_lp, "subtour");
  if(rval) goto CLEANUP;
	  

  /* Build a row for each degree equation */
  for(int i = 0; i < m_graph.node_count; i++) {
    rval = PSEPlp_new_row (&m_lp, 'E', 2.0);
    if(rval) goto CLEANUP;
  }

  /* Build a column for each edge of the Graph */
  for(int j = 0; j < m_graph.edge_count; j++) {
    int *nodes = (int*)m_graph.edges[j].end;
    double objective_val = (double)m_graph.edges[j].len;
    rval = PSEPlp_addcols (&m_lp, num_vars, num_non_zero, &objective_val,
			   &cmatbeg, nodes, coefficients, &lower_bound,
			   &upper_bound);
    if(rval) goto CLEANUP;
  }

  prefs = _prefs;

  cout << "Adding at most " << prefs.max_per_round << " cuts per round, "
       << "keeping at most " << prefs.q_max_size << " candidate blossoms\n";


  try {
    m_lp_edges.resize(m_graph.edge_count);
    old_colstat.resize(m_graph.edge_count, CPX_AT_LOWER);
    old_rowstat.resize(m_graph.node_count, CPX_AT_LOWER);
    frac_colstat.resize(m_graph.edge_count);
    frac_rowstat.resize(m_graph.edge_count);
  } catch(const std::bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Problem allocating LP vectors, ");
  }

 CLEANUP:
  if(rval){
    cerr << "LPGroup constructor failed\n";
    m_lp_edges.clear();
    old_colstat.clear();
    old_rowstat.clear();
    frac_colstat.clear();
    frac_rowstat.clear();
  }
}
