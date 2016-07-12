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

  cout << "BRANCH VERSION: GETOPT!!!!" << endl;

  if(initial_parse(argc, argv, graph, tour_node_indices)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }

  TSP_Solver solver(graph, tour_node_indices);
  solver.basis_init();

  double start = PSEP_zeit();
  int stat;
  int num_seg, num_2match;
  int segval, matchval;
  int rounds = 0;
  int total_cuts = 0;

  cout << "Pivoting until optimality or no more cuts" << endl;

  while(true){
    rounds++;
    if(solver.pivot_until_change(&stat))
      break;

    if(stat == 3)
      break;
    if(stat == 2){
      if(solver.update_best_tour())
	break;
      else{
	cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "!!!AUGMENTED TOUR!!!!" << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!!" << endl;
	continue;
      }
    }

    if(solver.pivot_back())
      break;

    segval = solver.cutcall.segment(&num_seg);
    if(segval == 1)
      break;
    total_cuts += num_seg;

    matchval = solver.cutcall.blossom(250 - num_seg, &num_2match);
    if(matchval == 1)
      break;

    total_cuts += num_2match;

    cout << "Added " << num_seg << " segment cuts and "
	 << num_2match << " blossom inequalities" << endl;

    if(segval + matchval == 4)
      break;
  }

  if(stat != 3)
    cout << "Terminated due to lack of cutting planes after "
	 << rounds << " rounds of separation" << endl;
  cout << total_cuts << " cutting planes added over "
       << rounds << " rounds of separation" << endl;


  cout << "Finished with runtime " << PSEP_zeit() - start << endl;

  return 0;
}

static int initial_parse(int ac, char **av, Graph &graph,
			 vector<int> &node_indices){
  char *fname = (char *) NULL;
  int seed = 0;
  bool devex = false;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "ads:")) != EOF) {
    switch(c) {
    case 'd':
      devex = true;
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

  LP::devex_switch = devex;
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
  fprintf(stderr, "   -d    enable dynamic switch to devex\n");
  fprintf(stderr, "   -s d  set random seed to d\n");
}
