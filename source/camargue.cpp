#include "config.hpp"

#ifndef CMR_DO_TESTS

#include "solver.hpp"
#include "util.hpp"
#include "io_util.hpp"
#include "timer.hpp"

#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

using std::string;
using std::cout;
using std::cerr;

using std::unique_ptr;

using std::runtime_error;
using std::logic_error;
using std::exception;


static void initial_parse(int argc, char **argv,
                          string &tsp_fname, string &tour_fname,
                          int &seed, int &rnodes, int &rgrid,
                          CMR::OutPrefs &outprefs, bool &sparseflag,
                          bool &branchflag);

static void usage(const std::string &fname);

int main(int argc, char** argv) try
{
    string tsp_fname;
    string tour_fname;
    
    int seed = 0;

    int rand_nodes = 0;
    int rand_grid = 0;

    bool sparse = false;
    bool branch = true;
    
    CMR::OutPrefs outprefs;

    unique_ptr<CMR::Solver> tsp_solver;

    initial_parse(argc, argv, tsp_fname, tour_fname,
                  seed, rand_nodes, rand_grid, outprefs, sparse, branch);

    if (!tsp_fname.empty()) {
        if (tour_fname.empty())
            tsp_solver = CMR::util::make_unique<CMR::Solver>(tsp_fname, seed,
                                                             outprefs);
        else
            tsp_solver = CMR::util::make_unique<CMR::Solver>(tsp_fname,
                                                             tour_fname, seed,
                                                             outprefs);
    } else
        tsp_solver = CMR::util::make_unique<CMR::Solver>(seed, rand_nodes,
                                                         rand_grid, outprefs);

    tsp_solver->cut_sel.safeGMI = sparse;

    CMR::Timer t;
    t.start();

    if (branch)
        tsp_solver->abc(!sparse);
    else
        tsp_solver->cutting_loop(!sparse, true, true);

    t.stop();
    t.report(true);

    return 0;
    
} catch (const exception &e) {
    cerr << "Exception in Camargue main: " << e.what() << "\n";
    return 1;
}

static void initial_parse(int ac, char **av,
                          string &tsp_fname, string &tour_fname,
                          int &seed, int &rnodes, int &rgrid,
                          CMR::OutPrefs &outprefs, bool &sparseflag,
                          bool &branchflag)
{
    bool randflag = false;
    
    rgrid = 1000000;

    int c;

    if (ac == 1) {
        usage(av[0]);
        throw logic_error("No arguments specified");
    }

    while ((c = getopt(ac, av, "aPRSn:g:s:t:")) != EOF) {
        switch (c) {
        case 'P':
            branchflag = false;
            break;
        case 'R':
            randflag = true;
            break;
        case 'S':
            sparseflag = true;
            break;
        case 'n':
            rnodes = atoi(optarg);
            break;
        case 'g':
            rgrid = atoi(optarg);
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 't':
            tour_fname = optarg;
            break;
        case '?':
        default:
            usage(av[0]);
            throw logic_error("Erroneous argument");
        }
    }

    if (optind < ac)
        tsp_fname = av[optind++];

    if (optind != ac) {
        usage(av[0]);
        throw logic_error("Bad option count");
    }

    if (tsp_fname.empty() && rnodes <= 0) {
        usage(av[0]);
        throw logic_error("Must specify problem file or random nodecount");
    }

    if (randflag &&
        (!tsp_fname.empty() || !tour_fname.empty())) {
        usage(av[0]);
        throw logic_error("Cannot specify filenames and random prob");
    }

    if (tsp_fname.empty() && !tour_fname.empty()) {
        usage(av[0]);
        throw logic_error("Cannot specify tour without TSPLIB file.");
    }
}

static void usage(const std::string &fname)
{
    cerr << "Usage: " << fname << " [-see below-] [-prob_file-]\n";
    cerr << "\t\t FLAG OPTIONS\n"
         << "-P \t Pure primal cutting plane solution: do not branch.\n"
         << "-R \t Generate random problem.\n"
         << "   \t Notes:\t Incompatible with -t flag and TSPLIB file\n"
         << "   \t       \t Must be set to specify -n, -g below\n"
         << "-S \t Sparse mode: solve over a sparse edge set, no pricing\n"
         << "   \t Notes:\t Must be set for safe Gomory cuts to be used.\n\n";
    cerr << "\t\t PARAMETER OPTIONS (argument x)\n"
         << "-n \t random problem with x nodes\n"
         << "-g \t random problem gridsize x by x (1 million default)\n"
         << "-s \t random seed x used throughout code (current time default)\n"
         << "-t \t load starting tour from path x\n";
}

#endif //CMR_DO_TESTS
