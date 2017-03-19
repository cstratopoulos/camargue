#include "config.hpp"

#ifndef CMR_DO_TESTS

#include "solver.hpp"
#include "abc_nodesel.hpp"
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
using std::endl;

using std::unique_ptr;

using std::runtime_error;
using std::logic_error;
using std::exception;


static void initial_parse(int argc, char **argv,
                          string &tsp_fname, string &tour_fname,
                          int &seed, int &rnodes, int &rgrid,
                          CMR::OutPrefs &outprefs, bool &sparseflag,
                          bool &branchflag, int &node_sel);

static void usage(const std::string &fname);

int main(int argc, char** argv) try
{
    string tsp_fname = "";
    string tour_fname = "";

    int seed = 0;

    int rand_nodes = 0;
    int rand_grid = 0;

    int node_sel = 0;

    bool sparse = false;
    bool branch = true;

    CMR::OutPrefs outprefs;

    unique_ptr<CMR::Solver> tsp_solver;

    initial_parse(argc, argv, tsp_fname, tour_fname,
                  seed, rand_nodes, rand_grid, outprefs, sparse, branch,
                  node_sel);

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
    tsp_solver->cut_sel.localcuts = true;

    CMR::Timer t(tsp_solver->inst_info().problem_name() + " overall");
    t.start();

    bool do_price = !sparse;

    if (branch) {
        switch (node_sel) {
        case 0:
            cout << "Interleaved best-tour/best-bound search" << endl;
            tsp_solver->abc<CMR::ABC::InterBrancher>(do_price);
            break;
        case 1:
            cout << "Best-tour branching with LK tours" << endl;
            tsp_solver->abc<CMR::ABC::TourBrancher>(do_price);
            break;
        case 2:
            cout << "Best-bound with primal strong branch probes" << endl;
            tsp_solver->abc<CMR::ABC::BoundBrancher>(do_price);
            break;
        case 3:
            cout << "DFS branching" << endl;
            tsp_solver->abc<CMR::ABC::DFSbrancher>(do_price);
            break;
        default:
            throw logic_error("Unimplemented node selection rule");
        }
    } else {
        tsp_solver->cutting_loop(!sparse, true, true);
    }

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
                          bool &branchflag, int &node_sel)
{
    bool randflag = false;

    rgrid = 1000000;

    int c;

    if (ac == 1) {
        usage(av[0]);
        throw logic_error("No arguments specified");
    }

    while ((c = getopt(ac, av, "aEGPRSVXb:n:g:s:t:")) != EOF) {
        switch (c) {
        case 'E':
            outprefs.save_tour_edges = true;
            break;
        case 'G':
            outprefs.gif_tour = true;
            break;
        case 'P':
            branchflag = false;
            break;
        case 'R':
            randflag = true;
            break;
        case 'S':
            sparseflag = true;
            break;
        case 'V':
            outprefs.verbose = true;
            break;
        case 'X':
            outprefs.dump_xy = true;
            break;
        case 'b':
            node_sel = atoi(optarg);
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
            throw logic_error("Bad argument");
        }
    }

    if (optind < ac)
        tsp_fname = av[optind++];

    if (optind != ac) {
        usage(av[0]);
        throw logic_error("Bad option count");
    }

    if (node_sel < 0 || node_sel > 3)
        throw logic_error("Branching strategy (-b) must be 0, 1, 2, or 3");

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
         << "-E \t Write tour edges to file (in addition to nodes).\n"
         << "-G \t GIF output: write each new tour to a distinct file.\n"
         << "-P \t Pure primal cutting plane solution: do not branch.\n"
         << "-R \t Generate random problem.\n"
         << "   \t Notes:\t Incompatible with -t flag and TSPLIB file\n"
         << "   \t       \t Must be set to specify -n, -g below\n"
         << "-S \t Sparse mode: solve over a sparse edge set, no pricing\n"
         << "   \t Notes:\t Must be set for safe Gomory cuts to be used.\n"
         << "-V \t Verbose: print lots of messages.\n"
         << "-X \t Dump XY: If applicable, write x-y coords to file.\n"
         << "    \t Notes:\t Works for random instances and most TSPLIB "
         << "instances,\n"
         <<"\t\t except some with MATRIX or EXPLICIT coordinates.\n";
    cerr << "\n";
    cerr << "\t\t PARAMETER OPTIONS (argument x)\n"
         << "-b \t Branching strategy x (see below).\n"
         << "   \t 0\tInterleaved best-tour/best-bound every 10 nodes"
         << " (default).\n"
         << "   \t 1\tBest-tour using sparse LK tours.\n"
         << "   \t 2\tBest-bound using primal strong branch bounds.\n"
         << "   \t 3\tDepth-first search.\n"
         << "-n \t Random problem with x nodes\n"
         << "-g \t Random problem gridsize x by x (1 million default)\n"
         << "-s \t Random seed x used throughout code (current time default)\n"
         << "-t \t Load starting tour from path x" << endl;
}

#endif //CMR_DO_TESTS
