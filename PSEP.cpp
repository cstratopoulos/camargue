#include<iostream>
#include<string>
#include<iomanip>
#include<memory>
#include<vector>

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <cstring>

#include "tsp_solver.h"

using namespace std;

static int initial_parse(int ac, char **av, string &fname,
			 PSEP::RandProb &randprob, PSEP::LP::Prefs &prefs);
static void usage(const string &fname);

int main(int argc, char* argv[]){
  PSEP::LP::Prefs prefs;
  PSEP::RandProb randprob;
  unique_ptr<CCdatagroup> dat(new CCdatagroup);
  string probfile;

  if(initial_parse(argc, argv, probfile, randprob, prefs)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }

  double overall = PSEP::zeit();
  PSEP::TSPSolver solver(probfile, randprob, prefs, dat);
  dat.reset();
  
  if(solver.call(PSEP::SolutionProtocol::PURECUT))
    exit(1);
  cout << "                    everything: "
       << PSEP::zeit() - overall << "\n";
}

static int initial_parse(int ac, char **av, string &fname,
			 PSEP::RandProb &randprob,
			 PSEP::LP::Prefs &prefs){
  bool rand = false;
  int pricing_choice = 0;
  int dp_factor = -1;
  int seed = 0;
  int ncount = 0;
  int gridsize = 100;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "aD:p:Rg:n:s:")) != EOF) {
    switch(c) {
    case 'D':
      dp_factor = atoi(optarg);
      break;
    case 'p':
      pricing_choice = atoi(optarg);
      break;
    case 'R':
      rand = true;
      break;
    case 'g':
      gridsize = atoi(optarg);
      break;
    case 'n':
      ncount = atoi(optarg);
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

  if(optind < ac)
    fname = av[optind++];
  if(optind != ac){
    usage(av[0]);
    return 1;
  }

  if(fname.empty() && ncount == 0){
    printf("Must specify a problem file or nodecount for random prob\n");
    return 1;
  }

  if(!fname.empty() && rand){
    printf("Cannot specify both filename and random problem\n");
    return 1;
  }

  randprob.nodecount = ncount;
  randprob.gridsize = gridsize;
  randprob.seed = seed;

  switch(pricing_choice){
  case 0:
    prefs.price_method = PSEP::LP::Pricing::Devex;
    cout << "Devex pricing\n";
    break;
  case 1:
    prefs.price_method = PSEP::LP::Pricing::SlackSteepest;
    cout << "Steepest edge w slack initial norms\n";
    break;
  case 2:
    prefs.price_method = PSEP::LP::Pricing::Steepest;
    cout << "True steepest edge\n";
    break;
  default:
    cout << "Pricing method " << pricing_choice << " out of range\n";
    usage(av[0]);
    return 1;
  }

  prefs.dp_threshold = 5 * dp_factor;
  if(dp_factor >= 0)
    cout << "DP separation will be tried after "
	 << prefs.dp_threshold << " non-degenerate pivots w no augmentation\n";

  return 0;
}

static void usage(const string &fname){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", fname.data());
  fprintf(stderr, "-------FLAG OPTIONS ---------------------------------\n");
  fprintf(stderr, "-R    generate random problem\n");
  fprintf(stderr, "------ PARAMETER OPTIONS (argument x) ---------------\n");
  fprintf(stderr, "-D    only call simpleDP sep after 5x rounds of cuts \n");
  fprintf(stderr, "      with no augmentation. (disabled by default)\n");
  fprintf(stderr, "-p     set primal pricing protocol to:\n");
  fprintf(stderr, "    0 (default) devex\n");
  fprintf(stderr, "    1 steepest edge with slack initial norms\n");
  fprintf(stderr, "    2 true steepest edge.\n");
  fprintf(stderr, "-g    gridsize for random problem (100 default)\n");
  fprintf(stderr, "-n    nodecount for random problem (must be nonzero if\n");
  fprintf(stderr, "      -R is used\n");
  fprintf(stderr, "-s    random seed for random problem\n");
}
