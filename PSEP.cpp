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
static int initial_parse(int ac, char **av,char **fname, PSEP_LP_Prefs &prefs);
static void usage(char *f);

int main(int argc, char* argv[]){
  PSEP_LP_Prefs prefs;
  CCdatagroup *dat = new CCdatagroup;
  char *fname;

  cout << "BRANCH VERSION: MASTER (ABC in progress)\n";

  if(initial_parse(argc, argv, &fname, prefs)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }


  TSP_Solver solver(fname, prefs, dat);

  delete dat;

  return solver.call(PSEP::SolutionProtocol::PURECUT);
}

static int initial_parse(int ac, char **av, char **fname,
			 PSEP_LP_Prefs & prefs){
  *fname = (char *) NULL;
  int seed = 0;
  int pricing_choice = 0;
  int switching_choice = 0;
  int dp_factor = -1;
  bool jumpstart = false;
  bool redfix = true;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "ThxaD:d:jp:s:")) != EOF) {
    switch(c) {
    case 'T':
      tooth = 1;
      break;
    case 'x':
      redfix = false;
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

  if(!redfix){
    prefs.redcost_fixing = false;
    cout << "Edges will not be fixed via reduced cost, ";
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

  prefs.dp_threshold = 5 * dp_factor;
  if(dp_factor >= 0)
    cout << "DP separation will be tried after "
	 << prefs.dp_threshold << " non-degenerate pivots w no augmentation\n";
  
  UTIL::seed = seed;

  return 0;
}

static void usage(char *f){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", f);
  fprintf(stderr, "------ FLAG OPTIONS (no argument)--------------------\n");
  fprintf(stderr, "-T       enable simpleDP/simple tooth testing.\n");
  fprintf(stderr, "-j       jumpstart slow sequences of pivots with\n");
  fprintf(stderr, "         a temporary switch to steepest edge.\n");
  fprintf(stderr, "-x       turn off reduced cost fixing (on by default).\n");
  fprintf(stderr, "------ PARAMETER OPTIONS (argument x) ---------------\n");
  fprintf(stderr, "-D    only call simpleDP sep after 5x rounds of cuts \n");
  fprintf(stderr, "      with no augmentation. (disabled by default)\n");
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
}
