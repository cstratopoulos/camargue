/**
 * @file
 * @brief The camargue main program.
 */
#include "config.hpp"

#ifndef CMR_DO_TESTS

#include "solver.hpp"
#include "abc_nodesel.hpp"
#include "util.hpp"
#include "io_util.hpp"
#include "timer.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
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

static constexpr double large_neg{-std::numeric_limits<double>::max() + 1.0};

/// Data grabbed from the user at the command line.
struct OptData {
    string tsp_fname = "";
    string tour_fname = "";

    int seed = 0;

    int rand_nodes = 0;
    int rand_grid = 1000000;

    int node_sel = 0;
    int cut_sel = 1;
    int edge_sel = 0;

    bool sparse = false;
    bool branch = true;

    double target_lb{large_neg};
};


static void proc_label();
static void initial_parse(int argc, char **argv, OptData &opt_dat,
                          CMR::OutPrefs &outprefs);

static void usage(const std::string &fname);

int main(int argc, char** argv) try
{
    OptData opt_dat;
    CMR::OutPrefs outprefs;

    unique_ptr<CMR::Solver> tsp_solver;

    initial_parse(argc, argv, opt_dat, outprefs);

    using EdgePlan = CMR::Graph::EdgePlan;
    EdgePlan ep = ((opt_dat.edge_sel == 0) ? EdgePlan::Linkern :
                   EdgePlan::Delaunay);

    string &tsp_fname = opt_dat.tsp_fname;
    string &tour_fname = opt_dat.tour_fname;
    int &seed = opt_dat.seed;
    int &rand_nodes = opt_dat.rand_nodes;
    int &rand_grid = opt_dat.rand_grid;

    proc_label();

    if (!tsp_fname.empty()) {
        if (tour_fname.empty())
            tsp_solver = CMR::util::make_unique<CMR::Solver>(tsp_fname, seed,
                                                             ep,
                                                             outprefs);
        else
            tsp_solver = CMR::util::make_unique<CMR::Solver>(tsp_fname,
                                                             tour_fname, seed,
                                                             ep,
                                                             outprefs);
    } else
        tsp_solver = CMR::util::make_unique<CMR::Solver>(seed, rand_nodes,
                                                         rand_grid, ep,
                                                         outprefs);

    if (opt_dat.target_lb != large_neg)
        tsp_solver->set_lowerbound(opt_dat.target_lb);

    using SelPreset = CMR::Solver::CutSel::Presets;
    int &cut_sel = opt_dat.cut_sel;
    bool &sparse = opt_dat.sparse;

    if (cut_sel == 1) {
            tsp_solver->choose_cuts(sparse ? SelPreset::Sparse :
                                    SelPreset::Aggressive);
    } else if (cut_sel == 0) {
        tsp_solver->choose_cuts(SelPreset::Vanilla);
        tsp_solver->cut_sel.safeGMI = sparse;
    } else
        throw logic_error("Unimplemented CutSel preset val");

    CMR::Timer t(tsp_solver->inst_info().problem_name() + " overall");
    t.start();

    bool do_price = !sparse;

    if (opt_dat.branch) {
        switch (opt_dat.node_sel) {
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

static void initial_parse(int ac, char **av, OptData &opt_dat,
                          CMR::OutPrefs &outprefs)
{
    bool randflag = false;

    int c;

    if (ac == 1) {
        usage(av[0]);
        throw logic_error("No arguments specified");
    }

    while ((c = getopt(ac, av, "aBEGPRSVXb:c:e:l:n:g:s:t:")) != EOF) {
        switch (c) {
        case 'B':
            outprefs.prog_bar = true;
            break;
        case 'E':
            outprefs.save_tour_edges = true;
            break;
        case 'G':
            outprefs.gif_tour = true;
            break;
        case 'P':
            opt_dat.branch = false;
            break;
        case 'R':
            randflag = true;
            break;
        case 'S':
            opt_dat.sparse = true;
            break;
        case 'V':
            outprefs.verbose = true;
            break;
        case 'X':
            outprefs.dump_xy = true;
            break;
        case 'b':
            opt_dat.node_sel = atoi(optarg);
            break;
        case 'c':
            opt_dat.cut_sel = atoi(optarg);
            break;
        case 'e':
            opt_dat.edge_sel = atoi(optarg);
            break;
        case 'l':
            opt_dat.target_lb = atof(optarg);
            break;
        case 'n':
            opt_dat.rand_nodes = atoi(optarg);
            break;
        case 'g':
            opt_dat.rand_grid = atoi(optarg);
            break;
        case 's':
            opt_dat.seed = atoi(optarg);
            break;
        case 't':
            opt_dat.tour_fname = optarg;
            break;
        case '?':
        default:
            usage(av[0]);
            throw logic_error("Bad argument");
        }
    }

    if (optind < ac)
        opt_dat.tsp_fname = av[optind++];

    if (optind != ac) {
        usage(av[0]);
        throw logic_error("Bad option count");
    }

    if (opt_dat.node_sel < 0 || opt_dat.node_sel > 3) {
        usage(av[0]);
        throw logic_error("Branching strategy (-b) must be 0, 1, 2, or 3");
    }

    if (opt_dat.cut_sel < 0 || opt_dat.cut_sel > 1) {
        usage(av[0]);
        throw logic_error("Cut sel (-c) must be 0 or 1");
    }

    if (opt_dat.edge_sel < 0 || opt_dat.edge_sel > 1) {
        usage(av[0]);
        throw logic_error("Edge sel (-e) must be 0 or 1");
    }

    if (opt_dat.tsp_fname.empty() && opt_dat.rand_nodes <= 0) {
        usage(av[0]);
        throw logic_error("Must specify problem file or random nodecount");
    }

    if (randflag &&
        (!opt_dat.tsp_fname.empty() || !opt_dat.tour_fname.empty())) {
        usage(av[0]);
        throw logic_error("Cannot specify filenames and random prob");
    }

    if (opt_dat.tsp_fname.empty() && !opt_dat.tour_fname.empty()) {
        usage(av[0]);
        throw logic_error("Cannot specify tour without TSPLIB file.");
    }

    if (outprefs.verbose && outprefs.prog_bar) {
        usage(av[0]);
        throw logic_error("Requested progress bar and verbose.");
    }
}

static void proc_label()
{
    char buf[1024];
    ::gethostname(buf, 1024);
    printf("Hostname: %s Current PID: %d", buf, static_cast<int>(::getpid()));
    cout << endl;
}

static void usage(const std::string &fname)
{
    cerr << "Usage: " << fname << " [-see below-] [-prob_file-]\n";
    cerr << "\t\t FLAG OPTIONS\n"
         << "-B \t Show an 80-column progress bar for piv values "
         << "(incompatible w verbose).\n"
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
         << "-c \t Cut selection x (see below, each contains prev).\n"
         << "   \t 0\tPrimal algorithms + simple standard heuristics.\n"
         << "   \t 1\tAs above, plus local cuts/cut metamorphoses (default).\n"
         << "-e \t Initial edge set (see below).\n"
         << "   \t 0\tUnion of 10 LK tours (+ quad-2 on tiny probs) (default).\n"
         << "   \t 1\tEuclidean-norm Delaunay triangulation.\n"
         << "   \t Notes:\t If a Delaunay triangulation is requested with an\n"
         << "   \t incompatible norm, the Linkern edges will be used.\n"
         << "-g \t Random problem gridsize x by x (1 million default)\n"
         << "-l \t Target lower bound: report optimal if tour is at most x.\n"
         << "-n \t Random problem with x nodes\n"
         << "-s \t Random seed x used throughout code (current time default)\n"
         << "-t \t Load starting tour from path x" << endl;
}

#endif //CMR_DO_TESTS
