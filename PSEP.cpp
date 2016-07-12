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
			 vector<int> &node_indices);
static void usage(char *f);

int main(int argc, char* argv[]){
  Graph graph;
  vector<int> tour_node_indices;

  cout << "BRANCH VERSION: PRICING OPTIONS" << endl;

  if(initial_parse(argc, argv, graph, tour_node_indices)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }

  TSP_Solver solver(graph, tour_node_indices);

  double start = PSEP_zeit();

  solver.pure_cut();

  cout << "Finished with runtime " << PSEP_zeit() - start << endl;

  return 0;
}

static int initial_parse(int ac, char **av, Graph &graph,
			 vector<int> &node_indices){
  char *fname = (char *) NULL;
  int seed = 0;
  int dynamic_switch = 0;
  int pricing_method = 0;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "ad:p:s:")) != EOF) {
    switch(c) {
    case 'd':
      dynamic_switch = atoi(optarg);
      break;
    case 'p':
      pricing_method = atoi(optarg);
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

  if(dynamic_switch < 0 || dynamic_switch > 2){
    printf("Dynamic switch %d is out of range\n", dynamic_switch);
    usage(av[0]);
    return 1;
  }

  if(pricing_method < 0 || pricing_method > 2){
    printf("Pricing method %d is out of range\n", pricing_method);
    usage(av[0]);
    return 1;
  }

  LP::PRICING::SWITCHING::choice = dynamic_switch;
  LP::PRICING::choice = pricing_method;

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
  fprintf(stderr, "      x = 0 do not switch pricing methods\n");
  fprintf(stderr, "        = 1 switch when a non-degenerate pivot takes");
  fprintf(stderr, "            more than 3 * number of rows iterations\n");
  fprintf(stderr, "        = 2 switch from the start\n");
  fprintf(stderr, "   -p x   set primal pricing method to x\n");
  fprintf(stderr, "      x = 0 devex (default)\n");
  fprintf(stderr, "      x = 1 steepest edge, slack init norms\n");
  fprintf(stderr, "      x = 3 full-blown steepest edge\n");
  fprintf(stderr, "   -s x   set random seed to x\n");
}
