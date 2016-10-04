#include<iostream>
#include<string>
#include<iomanip>
#include<memory>
#include<vector>

#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<getopt.h>

#include "tsp_solver.hpp"
#include "setbank.hpp"
#include "PSEP_util.hpp"


using namespace std;

static int initial_parse(int ac, char **av, string &fname,
			 PSEP::RandProb &randprob, PSEP::LP::Prefs &prefs,
			 bool &sparseflag);

static void usage(const string &fname);

int main(int argc, char* argv[]){  
  PSEP::LP::Prefs prefs;
  PSEP::RandProb randprob;
  unique_ptr<CCdatagroup> dat(new CCdatagroup);
  string probfile;
  bool do_sparse = false;

  if(initial_parse(argc, argv, probfile, randprob, prefs, do_sparse)){
    cerr << "Problem parsing arguments" << endl;
    exit(1);
  }

  double overall = PSEP::zeit();
  PSEP::TSPSolver solver(probfile, randprob, prefs, dat, do_sparse);
  dat.reset();
  
  if(solver.call(PSEP::SolutionProtocol::PURECUT, do_sparse))
    exit(1);
  cout << "                    everything: "
       << PSEP::zeit() - overall << "\n";
}

static int initial_parse(int ac, char **av, string &fname,
			 PSEP::RandProb &randprob,
			 PSEP::LP::Prefs &prefs, bool &sparseflag){
  bool rand = false;
  int pricing_choice = 0;
  int dp_factor = -1;
  int cuts_per_round = 1;
  int max_q_size = 150;
  int seed = 0;
  int ncount = 0;
  int gridsize = 100;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "aD:c:q:p:RSg:n:s:")) != EOF) {
    switch(c) {
    case 'D':
      dp_factor = atoi(optarg);
      break;
    case 'c':
      cuts_per_round = atoi(optarg);
      break;
    case 'q':
      max_q_size = atoi(optarg);
      break;
    case 'p':
      pricing_choice = atoi(optarg);
      break;
    case 'R':
      rand = true;
      break;
    case 'S':
      sparseflag = true;
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
    cout << "DP separation will be tried every "
	 << prefs.dp_threshold << " non-degenerate pivots w no augmentation\n";

  if(cuts_per_round < 1 || max_q_size < 1 || max_q_size < cuts_per_round){
    cerr << "Invalid cuts per round or queue capacity\n";
    usage(av[0]);
    return 1;
  }
  prefs.max_per_round = cuts_per_round;
  prefs.q_max_size = max_q_size;

  return 0;
}
//TODO update gomory sep to conform to cuts per round
static void usage(const string &fname){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", fname.data());
  fprintf(stderr, "-------FLAG OPTIONS ------------------------------------\n");
  fprintf(stderr, "-R    generate random problem\n");
  fprintf(stderr, "-S    only solve sparse instance with edge set from \n");
  fprintf(stderr,"       10 Lin-Kernighan tours\n");
  fprintf(stderr, "------ PARAMETER OPTIONS (argument x) ------------------\n");
  fprintf(stderr, "-c    add at most x cuts per round (default 2)\n");
  fprintf(stderr, "-q    keep a queue of at most x blossom cuts to \n");
  fprintf(stderr, "      be checked before calling the blossom separation \n");
  fprintf(stderr, "      all over again (default 15)\n");
  fprintf(stderr, "-D    only call simpleDP sep after 5x rounds of cuts \n");
  fprintf(stderr, "      with no augmentation. (disabled by default)\n");
  fprintf(stderr, "-p    set primal pricing protocol to:\n");
  fprintf(stderr, "    0 (default) devex\n");
  fprintf(stderr, "    1 steepest edge with slack initial norms\n");
  fprintf(stderr, "    2 true steepest edge.\n");
  fprintf(stderr, "-g    gridsize for random problem (100 default)\n");
  fprintf(stderr, "-n    nodecount for random problem (must be nonzero if\n");
  fprintf(stderr, "      -R is used\n");
  fprintf(stderr, "-s    random seed for random problem and Lin-Kernighan.\n");
  fprintf(stderr, "      If not set, current time will be used.\n");
  fprintf(stderr, "      If set, initial tour will be constructed from only\n");
  fprintf(stderr, "      one run of Lin-Kernighan (rather than independent \n");
  fprintf(stderr, "      trials), to allow reproducibility.\n");
}
