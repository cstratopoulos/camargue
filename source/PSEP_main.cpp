/** \mainpage
 * This is the documentation mainpage for PSEPtsp, a TSP solver based on 
 * primal cutting plane methods. 
 *
 *
 * See @ref TSPSolver for information on the main solver class.
 *
 * See the Data namespace for information on internal data structures 
 * maintained and updated in the solution process
 *
 * See the  LP namespace for classes, constants, and structs related to 
 * managing the %LP relaxations in the solution process.
 *
 * See PureCut or ABC for overviews of the two main solution protocols. 
 *
 */

#include "tsp_solver.hpp"
#include "PSEP_util.hpp"
#include "graph_io.hpp"

#include <iostream>
#include <string>
#include <iomanip>
#include <memory>
#include <utility>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

static int initial_parse(int ac, char **av, std::string &fname,
			 std::string &tour_fname,
			 PSEP::RandProb &randprob, PSEP::LP::Prefs &prefs,
			 PSEP::OutPrefs &o_prefs,
			 bool &sparseflag, int &qnearest);

static void usage(const std::string &fname);

int main(int argc, char* argv[]){
  PSEP::LP::Prefs prefs;
  PSEP::OutPrefs o_prefs;
  PSEP::RandProb randprob;
  std::unique_ptr<CCdatagroup> dat(new CCdatagroup);
  /**@todo make this a regular ptr bc it confuses valgrind maybe */
  std::string probfile, tourfile;
  bool do_sparse = false;
  int qnearest = 0;
  /**@todo probably put this somewhere else */

  if(initial_parse(argc, argv, probfile, tourfile, randprob, prefs, o_prefs,
  		   do_sparse, qnearest)){
    std::cerr << "Problem parsing arguments\n";
    exit(1);
  }

  double overall = PSEP::zeit();
  std::unique_ptr<PSEP::TSPSolver> solver;

  try {
    if(tourfile.empty())
      solver = PSEP::make_unique<PSEP::TSPSolver>(probfile, randprob,
						      o_prefs, prefs, dat,
						      do_sparse, qnearest);
    else
      solver = PSEP::make_unique<PSEP::TSPSolver>(probfile, tourfile,
						  o_prefs, prefs, dat,
						  do_sparse, qnearest);
  } catch(...) {
    return 1;
  }
      
  int rval = solver->call(PSEP::SolutionProtocol::PURECUT, do_sparse);

  std::cout << "                    everything: "
	    << PSEP::zeit() - overall << "\n";
  return rval;
}

static int initial_parse(int ac, char **av, std::string &fname,
			 std::string &tourfname,
			 PSEP::RandProb &randprob,
			 PSEP::LP::Prefs &prefs, PSEP::OutPrefs &o_prefs,
			 bool &sparseflag,
			 int &qnearest){
  bool rand = false;
  int pricing_choice = 0;
  int dp_factor = -1;
  int cuts_per_round = 4;
  int max_q_size = 150;
  int seed = 0;
  int ncount = 0;
  int gridsize = 100;
  int tourprefs = 0;
  qnearest = 0;

  int c;

  if(ac == 1){
    usage(av[0]);
    return 1;
  }

  while((c = getopt(ac, av, "ad:c:q:p:RSXg:n:s:t:u:o:")) != EOF) {
    switch(c) {
    case 'd':
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
    case 'X':
      o_prefs.dump_xy = true;
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
    case 'u':
      qnearest = atoi(optarg);
      break;
    case 'o':
      tourprefs = atoi(optarg);
      break;
    case 't':
      tourfname = optarg;
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
    usage(av[0]);
    return 1;
  }

  if(!fname.empty() && rand){
    printf("Cannot specify both filename and random problem\n");
    usage(av[0]);
    return 1;
  }

  randprob.nodecount = ncount;
  randprob.gridsize = gridsize;
  randprob.seed = seed;

  switch(pricing_choice){
  case 0:
    prefs.price_method = PSEP::LP::Pricing::Devex;
    std::cout << "Devex pricing\n";
    break;
  case 1:
    prefs.price_method = PSEP::LP::Pricing::SlackSteepest;
    std::cout << "Steepest edge w slack initial norms\n";
    break;
  case 2:
    prefs.price_method = PSEP::LP::Pricing::Steepest;
    std::cout << "True steepest edge\n";
    break;
  default:
    std::cout << "Pricing method " << pricing_choice << " out of range\n";
    usage(av[0]);
    return 1;
  }

  prefs.dp_threshold = 5 * dp_factor;
  if(dp_factor >= 0)
    std::cout << "DP separation will be tried every "
	 << prefs.dp_threshold << " non-degenerate pivots w no augmentation\n";

  if(cuts_per_round < 1 || max_q_size < 1 || max_q_size < cuts_per_round){
    std::cerr << "Invalid cuts per round or queue capacity\n";
    usage(av[0]);
    return 1;
  }
  prefs.max_per_round = cuts_per_round;
  prefs.q_max_size = max_q_size;

  if(qnearest < 0 || qnearest > 10){
    std::cerr << "Invalid choice of quad-nearest density\n";
    usage(av[0]);
    return 1;
  }

  switch(tourprefs) {
  case 3:
    o_prefs.save_tour = false;
    o_prefs.save_tour_edges = false;
    break;
  case 1:
    o_prefs.save_tour = false;
    o_prefs.save_tour_edges = true;
    break;
  case 2:
    o_prefs.save_tour = true;
    o_prefs.save_tour_edges = true;
    break;
  case 0: default:
    o_prefs.save_tour = true;
    break;
  }

  return 0;
}

static void usage(const std::string &fname){
  fprintf(stderr, "Usage: %s [-see below-] [prob_file]\n", fname.data());
  fprintf(stderr, "-------FLAG OPTIONS ------------------------------------\n");
  fprintf(stderr, "-R    generate random problem\n");
  fprintf(stderr, "-S    only solve sparse instance with edge set from \n");
  fprintf(stderr,"       10 Lin-Kernighan tours\n");
  fprintf(stderr, "-X    dump xy-coords to file, if possible \n");
  fprintf(stderr, "------ PARAMETER OPTIONS (argument x) ------------------\n");
  fprintf(stderr, "-c    add at most x cuts per round (default 2)\n");
  fprintf(stderr, "-d    only call simpleDP sep after 5x rounds of cuts \n");
  fprintf(stderr, "      with no augmentation. (disabled by default)\n");
  fprintf(stderr, "-g    gridsize for random problem (100 default)\n");
  fprintf(stderr, "-n    nodecount for random problem (must be nonzero if\n");
  fprintf(stderr, "      -R is used\n");
  fprintf(stderr, "-o    tour file output level \n");
  fprintf(stderr, "    0 (default) only nodes of best tour\n");
  fprintf(stderr, "    1 only edges of best tour, node node format\n");
  fprintf(stderr, "    2 both nodes and edges of best tour\n");
  fprintf(stderr, "    3 nothing about best tour\n");
  fprintf(stderr, "-p    set primal pricing protocol to:\n");
  fprintf(stderr, "    0 (default) devex\n");
  fprintf(stderr, "    1 steepest edge with slack initial norms\n");
  fprintf(stderr, "    2 true steepest edge.\n");
  fprintf(stderr, "-q    keep a queue of at most x blossom cuts to \n");
  fprintf(stderr, "      be checked before calling the blossom separation \n");
  fprintf(stderr, "      all over again (default 15)\n");
  fprintf(stderr, "-s    random seed for random problem and Lin-Kernighan.\n");
  fprintf(stderr, "      If not set, current time will be used.\n");
  fprintf(stderr, "-t    load initial tour in filename x.\n");
  fprintf(stderr, "-u    initial edge set will be union of 10 LK tours plus\n");
  fprintf(stderr, "      quad x-nearest edges (0 default).\n");
}
