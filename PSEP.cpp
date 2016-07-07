#include<iostream>
#include<iomanip>
#include<vector>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <cstring>

#include "TSP_solver.h"

using namespace std;

void print_tour(const vector<int> &tour_indices, const Graph &graph);
static int load_tsplib (Graph &graph, CCdatagroup *dat, int ac, char **av);
static int initialize_lk_tour (Graph &graph, CCdatagroup *dat,
			       vector<int> &node_indices);

int main(int argc, char* argv[]){
  Graph graph;
  vector<int> tour_node_indices;
  CCdatagroup dat;

  cout << "BRANCH VERSION: Master" << endl;

  if(load_tsplib(graph, &dat, argc, argv)){
    cerr << "Problem getting tsplib" << endl;
    exit(1);
  }
  
  if(initialize_lk_tour(graph, &dat, tour_node_indices)){
    cerr << "Problem getting LK tour" << endl;
    exit(1);
  }

  TSP_Solver solver(graph, tour_node_indices);
  solver.print_best_tour_nodes();
  solver.basis_init();

  int stat = 1;
  int old_basic, old_nb, old_stat;

  cout << "Pivoting until no more segment cuts" << endl;
  int num_added = -1;

  double start = PSEP_zeit();

  while(true){
    if(solver.pivot_until_change(&old_basic, &old_nb, &old_stat, &stat))
      break;

    cout << "Pivot status: ";
    switch(stat){
    case(0):
      cout << "Fractional" << endl;
      break;
    case(1):
      cout << "Integral subtour" << endl;
      break;
    case(2):
      cout << "New tour" << endl;
      break;
    case(3):
      cout << "Tour fathomed optimal" << endl;
    }

    if(stat == 2 || stat == 3)
      break;
    
    if(solver.blossom_cutcall(80, &num_added) == 1)
      break;
    if(num_added == 0){
      cout << "No blossom cuts found, loop terminating" << endl;
      break;
    }
    if(solver.pivot_back(old_basic, old_nb, old_stat))
      break;    
  }

  if(stat == 2 || stat == 3)
    cout << "Pivoted to new tour, nothing to do" << endl;
  else {
    cout << "Insert: Pivot again and call blossom separation" << endl;
    /*
    solver.blossom_cutcall(80, &num_added);
    solver.pivot_back(old_basic, old_nb, old_stat);
    solver.pivot_until_change(&old_basic, &old_nb, &old_stat, &stat);

    cout << "Pivot status: ";
    switch(stat){
    case(0):
      cout << "Fractional" << endl;
      break;
    case(1):
      cout << "Integral subtour" << endl;
      break;
    case(2):
      cout << "New tour" << endl;
      break;
    case(3):
      cout << "Tour fathomed optimal" << endl;
    }
    */
  }


  cout << "Finished with total runtime: " << PSEP_zeit() - start << endl;

  return 0;
}

static int load_tsplib (Graph &graph, CCdatagroup *dat, int ac, char **av){
  char *fname = (char *) NULL;
  int rval = 0;

  if (ac != 2){
    fprintf(stderr, "Enter name of TSPLIB file to be read\n");
    return 1;
  } else {
    fname = av[1];
  }

  CCutil_init_datagroup(dat);
  rval = CCutil_gettsplib(fname, &(graph.node_count), dat);
  if (rval){
    fprintf(stderr, "get tsplib failed\n");
    return rval;
  }
  
  int ncount = graph.node_count;
  int ecount = (ncount * (ncount -1 )) / 2;
  graph.edge_count = ecount;
  graph.edges.resize(ecount);
  int e_index = 0;
  
  for (int i = 0; i < ncount; i++){
    for (int j = i+1; j < ncount; j++){
      graph.edges[e_index].end[0] = i;
      graph.edges[e_index].end[1] = j;
      graph.edges[e_index].len = CCutil_dat_edgelen(i, j, dat);
      e_index++;
    }
  }

  return rval;
}

static int initialize_lk_tour (Graph &graph, CCdatagroup *dat,
			       vector<int> &node_indices){
  int rval = 0;
  CCrandstate rand_state;
  CCedgegengroup plan;
  int ncount = graph.node_count;
  int ecount = 0;
  int *elist = (int *) NULL;
  int tcount = 0;
  int *tlist = (int *) NULL;
  //int *bestcyc = (int *) NULL;
  int *perm = (int *) NULL;
  int *cyc = (int *) NULL;
  double bestval, szeit;
  int silent = 1;
  int kicks = (ncount > 400 ? 100 : ncount / 4);
  int istour;
  int seed = (int) PSEP_real_zeit();
  seed = 1466719036; //breaks dsj1000 and pr1002
  //seed = 1466720112; //breaks pr2392
  cout << "LK seed: " << seed << endl;

  szeit = CCutil_zeit ();
  bestval = CCtsp_LP_MAXDOUBLE;

  //code copies from static int find_tour from concorde
  CCutil_sprand(seed, &rand_state);
  CCrandstate *rstate = &rand_state;
  perm = CC_SAFE_MALLOC (ncount, int);
  CCcheck_NULL (perm, "out of memory for perm");
  
  node_indices.resize(ncount);

  CCedgegen_init_edgegengroup (&plan);
  plan.quadnearest = 2;
  rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &ecount,
			  &elist, silent, rstate);
  CCcheck_rval (rval, "CCedgegen_edges failed");
  plan.quadnearest = 0;

  plan.tour.greedy = 1;
  rval = CCedgegen_edges (&plan, ncount, dat, (double *) NULL, &tcount,
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

  rval = CClinkern_tour (ncount, dat, ecount, elist, ncount, kicks,
			 cyc, &node_indices[0], &bestval, silent, 0.0, 0.0,
			 (char *) NULL,
			 CC_LK_GEOMETRIC_KICK, rstate);
  CCcheck_rval (rval, "CClinkern_tour failed");
  //end of copied code
  
  ecount = graph.edge_count;
  cout << "bestcyc has value " << bestval << endl;

 CLEANUP:

  CC_IFFREE (cyc, int);
  CC_IFFREE (elist, int);
  CC_IFFREE (tlist, int);
  return rval;
}

void print_tour(const vector<int> &tour_indices, const Graph &graph) {
    double tour_length = 0.0;
    cout << "Optimal tour:" << endl;
    for(int i = 0; i < (int)tour_indices.size(); i++) {
        Edge e = graph.edges[tour_indices[i]];
        tour_length += e.len;
        cout << e.end[0] << " " << e.end[1] << " " << endl;
    }
    cout << "Optimal tour value: " << tour_length << endl;
}
