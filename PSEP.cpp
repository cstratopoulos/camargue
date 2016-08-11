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

static int initial_parse(int ac, char **av,char **fname,
			 PSEP::LP::Prefs &prefs);
static void usage(char *f);

int main(int argc, char* argv[]){
  PSEP::LP::Prefs prefs;
  CCdatagroup *dat = new CCdatagroup;
  char *fname;

  cout << "BRANCH VERSION: MASTER\n";

  if(initial_parse(argc, argv, &fname, prefs)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }


  PSEP::TSPSolver solver(fname, prefs, dat);

  delete dat;

  return solver.call(PSEP::SolutionProtocol::PURECUT);
}

static int initial_parse(int ac, char **av, char **fname,
			 PSEP::LP::Prefs &prefs){
  *fname = (char *) NULL;
  int pricing_choice = 0;
  int dp_factor = -1;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "aD:p:")) != EOF) {
    switch(c) {
    case 'D':
      dp_factor = atoi(optarg);
      break;
    case 'p':
      pricing_choice = atoi(optarg);
      break;
    case '?':
    default:
      usage(av[0]);
      return 1;
    }
  }

  if(optind < ac)
    *fname = av[optind++];
  if(optind != ac){
    usage(av[0]);
    return 1;
  }

  if(!(*fname)){
    printf("Must specify a problem file\n");
    return 1;
  }

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

static void usage(char *f){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", f);
  fprintf(stderr, "------ PARAMETER OPTIONS (argument x) ---------------\n");
  fprintf(stderr, "-D    only call simpleDP sep after 5x rounds of cuts \n");
  fprintf(stderr, "      with no augmentation. (disabled by default)\n");
  fprintf(stderr, "-p     set primal pricing protocol to:\n");
  fprintf(stderr, "    0 (default) devex\n");
  fprintf(stderr, "    1 steepest edge with slack initial norms\n");
  fprintf(stderr, "    2 true steepest edge.\n");
}
