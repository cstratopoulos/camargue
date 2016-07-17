#include<iostream>
#include<iomanip>
#include<vector>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <cstring>

#include "TSP_solver.h"

using namespace std;

static int load_tsplib (Graph &graph, CCdatagroup *dat, char *fname);
static int initialize_lk_tour (Graph &graph, CCdatagroup *dat,
			       vector<int> &node_indices);
static int initial_parse(int ac, char **av, Graph &graph,
			 vector<int> &node_indices, PSEP_LP_Prefs &prefs);
static void usage(char *f);

int main(int argc, char* argv[]){
  Graph graph;
  PSEP_LP_Prefs prefs;
  vector<int> tour_node_indices;

  cout << "BRANCH VERSION: SHARED_PTR" << endl;

  if(initial_parse(argc, argv, graph, tour_node_indices, prefs)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }

  TSP_Solver solver(graph, tour_node_indices, prefs);

  double start = PSEP_zeit();

  //  solver.simple_test();
  solver.pure_cut();

  cout << "Finished with runtime " << PSEP_zeit() - start << endl;

  return 0;
}

static int initial_parse(int ac, char **av, Graph &graph,
			 vector<int> &node_indices, PSEP_LP_Prefs &prefs){
  char *fname = (char *) NULL;
  int seed = 0;
  int pricing_choice = 0;
  int switching_choice = 0;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "ad:p:s:")) != EOF) {
    switch(c) {
    case 'd':
      switching_choice = atoi(optarg);
      break;
    case 'p':
      pricing_choice = atoi(optarg);
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case '?':
    default:
      usage(av[0]);
      return 1;
    }
  }

  if(optind < ac) fname = av[optind++];

  if(optind != ac){
    usage (av[0]);
    return 1;
  }

  if(!fname){
    printf ("Must specify a problem file\n");
    return 1;
  }

  switch(switching_choice){
  case 0:
    prefs.switching_choice = LP::PRICING::SWITCHING::OFF;
    cout << "No switching" << endl;
    break;
  case 1:
    prefs.switching_choice = LP::PRICING::SWITCHING::DYNAMIC;
    cout << "Dynamic switching" << endl;
    break;
  case 2:
    prefs.switching_choice = LP::PRICING::SWITCHING::START;
    cout << "Switch from the start" << endl;
    break;
  default:
    cout << "Switching method " << switching_choice << " out of range\n";
    usage(av[0]);
    return 1;
  }

  switch(pricing_choice){
  case 0:
    prefs.pricing_choice = LP::PRICING::DEVEX;
    cout << "Devex pricing\n";
    break;
  case 1:
    prefs.pricing_choice = LP::PRICING::STEEPEST;
    cout << "Steepest edge w slack initial norms\n";
    break;
  case 2:
    prefs.pricing_choice = LP::PRICING::STEEPEST_REAL;
    cout << "True steepest edge\n";
    break;
  default:
    cout << "Pricing method " << pricing_choice << " out of range\n";
    usage(av[0]);
    return 1;
  }
  
  UTIL::seed = seed;

  CCdatagroup dat;
  if(load_tsplib(graph, &dat, fname))
    return 1;

  if(initialize_lk_tour(graph, &dat, node_indices))
    return 1;

  return 0;
}

static int load_tsplib (Graph &graph, CCdatagroup *dat, char *fname){
  int rval = 0;

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
  int seed;
  if(UTIL::seed)
    seed = UTIL::seed;
  else
    seed = (int) PSEP_real_zeit();
  
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

static void usage(char *f){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", f);
  fprintf(stderr, "   -d x   set dynamic pricing switch behavior to x\n");
  fprintf(stderr, "      0 = do not switch pricing methods\n");
  fprintf(stderr, "      1 = switch when a non-degenerate pivot takes\n");
  fprintf(stderr, "          more than 3 * number of nodes iterations\n");
  fprintf(stderr, "      2 = switch from the start\n");
  fprintf(stderr, "   -p x   set primal pricing method to x\n");
  fprintf(stderr, "      0 = devex (default)\n");
  fprintf(stderr, "      1 = steepest edge, slack init norms\n");
  fprintf(stderr, "      2 = full-blown steepest edge\n");
  fprintf(stderr, "   -s x   set random seed to x\n");
}
