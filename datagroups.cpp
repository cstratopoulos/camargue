extern "C" {
#include <concorde/INCLUDE/linkern.h>
#include <concorde/INCLUDE/util.h>
#include <concorde/INCLUDE/edgegen.h>
}

#include "datagroups.h"
#include "mip.h"

using namespace std;
using namespace PSEP::Data;

GraphGroup::GraphGroup(const string &fname, unique_ptr<CCdatagroup> &dat){
  int rval = 0;
  CCdatagroup *rawdat = dat.get();
  char *filestring = const_cast<char *>(fname.data());

  CCutil_init_datagroup(rawdat);
  rval = CCutil_gettsplib(filestring, &(m_graph.node_count), rawdat);
  if (rval){
    fprintf(stderr, "get tsplib failed\n");
    exit(1);
  }
  
  int ncount = m_graph.node_count;
  int ecount = (ncount * (ncount -1 )) / 2;
  m_graph.edge_count = ecount;
  m_graph.edges.resize(ecount);
  int e_index = 0;
  
  for (int i = 0; i < ncount; i++){
    for (int j = i+1; j < ncount; j++){
      m_graph.edges[e_index].end[0] = i;
      m_graph.edges[e_index].end[1] = j;
      m_graph.edges[e_index].len = CCutil_dat_edgelen(i, j, rawdat);
      m_graph.edge_lookup.emplace(IntPair(i,j), e_index);
      e_index++;
    }
  }

  island.resize(m_graph.node_count);
  delta.resize(m_graph.edge_count, 0);
  edge_marks.resize(m_graph.node_count, 0);
}

BestGroup::BestGroup(const Graph &m_graph, unique_ptr<CCdatagroup> &dat){
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
  double bestval, val, szeit;
  int trials = 2 * (ceil((int) ncount / 100));
  int silent = 1;
  int kicks = 5 * ncount;
  int istour;
  int seed;
  seed = (int) PSEP_real_zeit();
  
  cout << "LK seed: " << seed << ", " << trials << " trials\n";

  szeit = CCutil_zeit ();
  bestval = INFINITY;

  //code copies from static int find_tour from concorde
  CCutil_sprand(seed, &rand_state);
  CCrandstate *rstate = &rand_state;
  cyc = CC_SAFE_MALLOC(ncount, int); //commented out to allow dummy tour
  if(!cyc){
    cerr << "Out of memory for find_tour\n";
    rval = 1; goto CLEANUP;
  }
  
  best_tour_nodes.resize(ncount);
  perm.resize(ncount);
  best_tour_edges.resize(m_graph.edge_count, 0);

  CCedgegen_init_edgegengroup (&plan);
  plan.quadnearest = 2;
  rval = CCedgegen_edges (&plan, ncount, rawdat, (double *) NULL, &ecount,
			  &elist, silent, rstate);
  CCcheck_rval (rval, "CCedgegen_edges failed");
  plan.quadnearest = 0;

  plan.tour.greedy = 1;
  rval = CCedgegen_edges (&plan, ncount, rawdat, (double *) NULL, &tcount,
			  &tlist, silent, rstate);
  CCcheck_rval (rval, "CCedgegen_edges failed");

  if (tcount != ncount) {
    fprintf (stderr, "wrong edgeset from CCedgegen_edges\n");
    rval = 1; goto CLEANUP;
  }

  rval = CCutil_edge_to_cycle (ncount, tlist, &istour, cyc);
  CCcheck_rval (rval, "CCutil_edge_to_cycle failed");
  if (istour == 0) {
    fprintf (stderr, "Starting tour has an error\n");
    rval = 1; goto CLEANUP;
  }
  CC_FREE (tlist, int);

  rval = CClinkern_tour (ncount, rawdat, ecount, elist, ncount, kicks,
			 cyc, &best_tour_nodes[0], &bestval, silent, 0.0, 0.0,
			 (char *) NULL,
			 CC_LK_GEOMETRIC_KICK, rstate);
  CCcheck_rval (rval, "CClinkern_tour failed");
  if(rval) goto CLEANUP;
  //end of copied code (from find_tour)
  cout << "LK initial run: " << bestval << endl;

  //begin copied code from find_good_tour in tsp_call.c
  //must be commented out for dummy tour
  
  for(int i = 0; i < trials; i++){
    rval = CClinkern_tour(ncount, rawdat, ecount, elist, ncount, kicks,
			  (int *) NULL, cyc, &val, 1, 0.0, 0.0,
			  (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate);
    if(rval)
      cerr << "CClinkern_tour failed\n";
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
			  CC_LK_GEOMETRIC_KICK, rstate);
    if(rval){
      cerr << "CClinkern_tour failed\n"; goto CLEANUP;
    }

    cout << "LK run from best tour: " << bestval << "\n";
    for(int j = 0; j < ncount; j++)
      best_tour_nodes[j] = cyc[j];
  }

  for(int i = 0; i < m_graph.node_count; i++)
    perm[best_tour_nodes[i]] = i;
  
  ecount = m_graph.edge_count;
  min_tour_value = bestval;

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
    }
  }

 CLEANUP:
  CC_IFFREE (cyc, int);
  CC_IFFREE (elist, int);
  CC_IFFREE (tlist, int);
  if(rval)
    exit(1);
}

LPGroup::LPGroup(const Graph &m_graph, PSEP::LP::Prefs &_prefs,
			   const vector<int> &perm){
  //Build the basic LP
  PSEPlp_init (&m_lp);
  PSEPlp_create (&m_lp, "subtour");

  //Set MIP/presolve parameters
  PSEPmip_param(&m_lp);

  /* Build a row for each degree equation */
  for(int i = 0; i < m_graph.node_count; i++) {
    PSEPlp_new_row (&m_lp, 'E', 2.0);
  }

  /* Build a column for each edge of the Graph */
  int cmatbeg = 0, num_vars = 1, num_non_zero = 2;
  double coefficients[2] = {1.0, 1.0};
  double lower_bound = 0.0;
  double upper_bound = 1.0;
  for(int j = 0; j < m_graph.edge_count; j++) {
    int *nodes = (int*)m_graph.edges[j].end;
    double objective_val = (double)m_graph.edges[j].len;
    PSEPlp_addcols (&m_lp, num_vars, num_non_zero, &objective_val,
		    &cmatbeg, nodes, coefficients, &lower_bound,
		    &upper_bound);
  }

  m_lp_edges.resize(m_graph.edge_count);
  prefs = _prefs;

  old_colstat.resize(m_graph.edge_count, CPX_AT_LOWER);
  old_rowstat.resize(m_graph.node_count, CPX_AT_LOWER);
  frac_colstat.resize(m_graph.edge_count);
  frac_rowstat.resize(m_graph.edge_count);
}
