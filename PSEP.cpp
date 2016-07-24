#include<iostream>
#include<string>
#include<iomanip>
#include<vector>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <cstring>

#include "TSP_solver.h"

using namespace std;

static int tooth = 0;
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

  cout << "BRANCH VERSION: MASTER\n";

  if(initial_parse(argc, argv, graph, tour_node_indices, prefs)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }

  TSP_Solver solver(graph, tour_node_indices, prefs);

  double start = PSEP_zeit();

  if(tooth)
    solver.simple_test();
  else
    solver.pure_cut();

  cout << "Finished with runtime " << PSEP_zeit() - start << endl;

  return 0;
}

static int initial_parse(int ac, char **av, Graph &graph,
			 vector<int> &node_indices,
			 PSEP_LP_Prefs & prefs){
  char *fname = (char *) NULL;
  int seed = 0;
  int pricing_choice = 0;
  int switching_choice = 0;
  int dp_factor = 3;
  bool jumpstart = false;
  int max_per_round = 2;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

while((c = getopt(ac, av, "Tam:D:d:jp:s:")) != EOF) {
    switch(c) {
    case 'T':
      tooth = 1;
      break;
    case 'D':
      dp_factor = atoi(optarg);
      break;
    case 'd':
      switching_choice = atoi(optarg);
      break;
    case 'j':
      jumpstart = true;
      break;
    case 'p':
      pricing_choice = atoi(optarg);
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'm':
      max_per_round = atoi(optarg);
      break;
    case '?':
    default:
      usage(av[0]);
      return 1;
    }
  }

  if(optind < ac)
    fname = av[optind++];
  if(optind != ac){
    usage(av[0]);
    return 1;
  }

  if(!fname){
    printf("Must specify a problem file\n");
    return 1;
  }

  switch(switching_choice){
  case 0:
    prefs.switching_choice = LP::PRICING::SWITCHING::OFF;
    cout << "No primal pricing switching\n";
    break;
  case 1:
    prefs.switching_choice = LP::PRICING::SWITCHING::DYNAMIC;
    cout << "Dynamic switching, ";
    break;
  case 2:
    prefs.switching_choice = LP::PRICING::SWITCHING::START;
    cout << "Switch from the start, ";
    break;
  default:
    cout << "Switching method " << switching_choice << " out of range\n";
    usage(av[0]);
    return 1;
  }

  if(jumpstart){
    prefs.jumpstart = true;
    cout << "Slow pivots will be jumpstarted, ";
  }

  switch(pricing_choice){
  case 0:
    prefs.pricing_choice = LP::PRICING::DEVEX;
    if(switching_choice)
      cout << "Devex pricing\n";
    break;
  case 1:
    prefs.pricing_choice = LP::PRICING::STEEPEST;
    if(switching_choice)
      cout << "Steepest edge w slack initial norms\n";
    break;
  case 2:
    prefs.pricing_choice = LP::PRICING::STEEPEST_REAL;
    if(switching_choice)
      cout << "True steepest edge\n";
    break;
  default:
    cout << "Pricing method " << pricing_choice << " out of range\n";
    usage(av[0]);
    return 1;
  }

  if(max_per_round < 2){
    cout << "Entered too few number of cuts per round\n";
    usage(av[0]);
    return 1;
  }

  cout << "at most " << max_per_round << "cuts will be added per round\n";
  prefs.max_cuts_round = max_per_round;

  if(dp_factor > 4){
    cout << "DP factor " << dp_factor << " too high\n";
    usage(av[0]);
    return 1;
  }
  prefs.dp_threshold = 5 * dp_factor;
  cout << "DP separation will be tried after "
       << prefs.dp_threshold << " non-degenerate pivots w no augmentation\n";
  
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
      graph.edge_lookup.emplace(IntPair(i,j), e_index);
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
  double bestval, val, szeit;
  int trials = 1;
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
  //  cyc = CC_SAFE_MALLOC(ncount, int); //commented out to allow dummy tour
  perm = CC_SAFE_MALLOC (ncount, int);
  if(!perm){
    cerr << "Out of memory for find_tour\n";
    rval = 1; goto CLEANUP;
  }
  
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
  if(rval) goto CLEANUP;
  //end of copied code (from find_tour)
  cout << "LK initial run: " << bestval << endl;

  /*
  //begin copied code from find_good_tour in tsp_call.c
  //must be commented out for dummy tour
  
  for(int i = 0; i < trials; i++){
    rval = CClinkern_tour(ncount, dat, ecount, elist, ncount, kicks,
			  (int *) NULL, cyc, &val, 1, 0.0, 0.0,
			  (char *) NULL, CC_LK_GEOMETRIC_KICK, rstate);
    if(rval)
      cerr << "CClinkern_tour failed\n";
    cout << "LK run " << i << ": " << val << "\n";
    if(val < bestval){
      for(int j = 0; j < ncount; j++)
	node_indices[j] = cyc[j];
      bestval = val;
    }
  }

  if (trials > 0){
    rval = CClinkern_tour(ncount, dat, ecount, elist, ncount, 2 * kicks,
			  &node_indices[0], cyc, &bestval, 1, 0.0, 0.0,
			  (char *) NULL,
			  CC_LK_GEOMETRIC_KICK, rstate);
    if(rval){
      cerr << "CClinkern_tour failed\n"; goto CLEANUP;
    }

    cout << "LK run from best tour: " << bestval << "\n";
    for(int j = 0; j < ncount; j++)
      node_indices[j] = cyc[j];
  }
  */
  
  ecount = graph.edge_count;


 CLEANUP:

  CC_IFFREE (cyc, int);
  CC_IFFREE (elist, int);
  CC_IFFREE (tlist, int);
  return rval;
}

static void usage(char *f){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", f);
  fprintf(stderr, "------ FLAG OPTIONS (no argument)--------------------\n");
  fprintf(stderr, "-T       enable simpleDP/simple tooth testing.\n");
  fprintf(stderr, "-j       jumpstart slow sequences of pivots with\n");
  fprintf(stderr, "         a temporary switch to steepest edge.\n");
  fprintf(stderr, "------ PARAMETER OPTIONS (argument x) ---------------\n");
  fprintf(stderr, "-D    only call simpleDP sep after 5x rounds\n");
  fprintf(stderr, "      of cuts with no augmentation.\n");
  fprintf(stderr, "-d    set dynamic switch of pricing protocol:\n");
  fprintf(stderr, "    0 (default) stay on reduced cost pricing\n");
  fprintf(stderr, "    1 switch after 3*ncount iterations needed\n");
  fprintf(stderr, "      for a non-degenerate pivot\n");
  fprintf(stderr, "    2 switch immediately\n");
  fprintf(stderr, "-p     set primal pricing protocol to:\n");
  fprintf(stderr, "    0 (default) devex\n");
  fprintf(stderr, "    1 steepest edge with slack initial norms\n");
  fprintf(stderr, "    2 true steepest edge.\n");
  fprintf(stderr, "-s    set LK random seed to x.\n");
  fprintf(stderr, "-m    set max number of cuts added per round\n");
  fprintf(stderr, "              of pivoting.\n");
}
