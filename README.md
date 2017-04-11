Camargue
========

This is the README for Camargue, a TSP solver based on primal
cutting plane methods. Camargue tries to move from tour to tour, either
proving that a tour is optimal or finding a better one. Camargue was
developed in tandem with the research in my master's thesis, "Primal
Cutting Plane Methods for the Traveling Salesman Problem".

The most effective TSP solvers, such
as [Concorde](http://www.math.uwaterloo.ca/tsp/concorde), are based on
a _dual fractional_ approach, which moves from one lower bound to the
next. Camargue is not competitive with Concorde, although it gives
respectable performance on small- and medium-sized instances. Rather,
Camargue functions as a proof of concept for primal cutting plane
methods. Its development has been a rich ground of testing and
experimentation, showing how well-studied dual fractional methods can
be adapted to handle large, highly degenerate linear programming
problems in the primal case.

To explore the code, a good starting point would be the header and
documentation for CMR::Solver. As the name suggests, this object
manages the primal cutting plane solution process. As for the
branching machinery, you may want to look at derived classes of
CMR::ABC::BaseBrancher: these are used to implement node selection
rules which guide primal branch-and-cut searches. For an even broader
view of branching machinery, look in the namespace CMR::ABC.


I have tried to document the source in a way that keeps header files
compact and readable. Except for class/structure definitions, the
documentation in a `.hpp` is mostly terse one-liners, with detailed
coumentation of function parameters and behaviors in the `.cpp`
files.

This page contains information on installing and using Camargue. You
can also view the GitHub pages version at

https://cstratopoulos.github.io/camargue/index.html

The text is mostly the same, but you will be able to navigate inline
to the Camargue documentation, created with
[Doxygen](http://www.stack.nl/~dimitri/doxygen/).


See below for

- [Installation instructions](#install)
- [Usage info](#usage)


[Installation](#install)
------------

We now discuss how to install Camargue. Camargue makes very heavy use
of two external dependencies. These are

- The TSP solver/combinatorial optimization library
[Concorde](http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm).

- The linear programming solver [CPLEX
12](http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud)

Both must be installed before proceeding. A Concorde installation with
CPLEX is best, but it should work with QSOPT too. For help, see

- [here](http://www.leandro-coelho.com/installing-concorde-tsp-with-cplex-linux/)
for Leandro Coelho's guide to doing an install with CPLEX; and
- [here](https://qmha.wordpress.com/2015/08/20/installing-concorde-on-mac-os-x/)
for David S. Johnson's guide to installing on Mac OS X with QSOPT; and
- the [Concorde README](http://www.math.uwaterloo.ca/tsp/concorde/DOC/README.html).

After installing Concorde, go into the directory `camargue/externals`
and create a symlink to the `concorde` directory. That is, `concorde`
must point to the folder containing the files `TSP`, `INCLUDE`,
`CUT`, etc.

After installing CPLEX, we take care of the `Makefile`. A template is
provided in the `scripts` folder; move it into the Camargue directory
with

    cp scripts/Makefile.template Makefile

(For users on the University of Waterloo Linux servers, do

    cp scripts/Makefile.template.UWlinux Makefile

These are already specified with a path to CPLEX and a C++ compiler,
so you can skip to the steps below about running the install script.)

Open the `Makefile` and edit the definitions `CPXDIR` and
`CPX_LIB`. Details and examples are given in the `Makefile`. Also, if
necessary, you can change the `CC :=` definition to specify a C++
compiler. Camargue is a C++11 application, so the compiler must be
compliant with (most of) the C++11 standard. Camargue has been developed and
tested on various Linux and Mac machines, using `g++ >= 4.7` (which
does _not_ support the standard in full) and `clang++ >= 3.8.0`.

No further edits to the Makefile should be necessary. After that, you
can run

    ./cmr_install.py

from the `camargue/` directory to configure and install Camargue.
The install script uses flag arguments to configure the installation to
your preferences, and to edit some other files. Running it with no
arguments, or with the flag `-h` or `--help` will print some usage
info, but we will do a bit more discussion here.

Information on external dependencies (other than Concorde and CPLEX)
is [here](externals/extdeps.md). When running `cmr_install.py`, the
simplest options are the two catch-all (or catch-none) flags, `-F,
--full` and `-B, --bare`. Doing

    ./cmr_install.py --full

will attempt to configure Camargue for a `full` install which uses all
of the external dependencies. If you are connected to the internet,
this will download and extract any externals not already present (or
symlink'd) in the `externals/` directory. This is recommended for best
performance, and to observe all the features described in my
thesis. The `--bare` option is the complementary flag: it will use no
externals besides Concorde and CPLEX. For an in-between approach, you
can explicitly specify one or more of `catch, safemir,` or `omp` with
the `-W` flag.

This script, and the ones it invokes, will generate a `config.hpp`
header, performing diagnostics on the specified C++ compiler and the
presence or absence of external dependencies. You can double check
`config.hpp`, and the `Makefile` (which may be edited too if `-W omp`
or `--full` is used) to make sure everything looks right.

After performing all these steps, you should be able to compile and
run Camargue by running `make` from the main directory, creating the
`camargue` executable.

Additionally, if you like, you can run the unit
tests, benchmarks, and experiments that I used to develop
Camargue (and write my thesis!) by using the recipe `make test`. This
requires downloading the [Catch](https://github.com/philsquared/Catch)
unit testing header is downloaded, and that you requested Catch in the
install script with either `--full` or `-W catch`. Install info is
[here](externals/extdeps.md), with specific
information on invoking the unit tests [here](source/tests/unittests.md).

[Usage](#usage)
------

This heading is about standard command line usage of Camargue. For
information on running tests/benchmarks,
see [their documentation](source/tests/unittests.md).

This section will try to elaborate a bit on the terse documentation
that you get from typing `./camargue` with no arguments.

Camargue accepts problems in two formats: TSPLIB instances with a
`.tsp` suffix, and randomly generated Euclidean instances. If you have
a folder called `problems` in the Camargue directory with TSPLIB
instances in it, you can run one with

    ./camargue problems/dantzig42.tsp

This will attempt to augment or prove the optimality of a starting
tour computed by Concorde's implementation of chained
Lin-Kernighan. You can also specify a starting tour with the '-t' flag:

    ./camargue problems/dantzig42.tsp -t test_data/tours/dantzig42.sol

will run the solver with the starting tour `dantzig42.sol`. The format
of solution files supported is _not_ the TSPLIB `.tour` format. Rather,
it should be a file whose first line is the instance node count, with
the following lines (with arbitrary spacing/indentation) giving
zero-indexed ordering of the nodes. If you would like to use a tour in the
TSPLIB `.tour` format, I have included a simple script for this
purpose. Given a tour file like `pr2392.tour`, running

    scripts/tour_to_sol.py pr2392.tour

will generate a file called `pr2392.sol` in the format just described.

When loading a starting tour, Camargue will make sure no obvious
mistakes are present, checking that no node appears twice and that the
indices are drawn from the proper range. You are free to specify an
abysmal starting tour, though!

Random problems can be generated with the flag argument `-R`, and some
additional arguments. To generate a 500-node instance on the 1,000 by
1,000 square grid, run

    ./camargue -Rn500 -g1000

So `-R` is the flag, `-n` specifies node count, and `-g` specifies the
gridsize.

For both styles of problems, Camargue will generate an initial edge
set consisting of the nodes in the tour found, plus the union of 10
quick Lin-Kernighan runs as implemented by Concorde's edge generation
code. For Euclidean-norm instances, the option `-e 1` can be used to
set the Delaunay triangulation edges as the starting edge set too.

Also for both styles of problems, you can pass a random seed with
`-s`. This is to allow reproducibility through all areas of the
code. (Note, however, that if OMP is [enabled](externals/extdeps.md),
non-determinism will still be present.)
For a random problem, this will be used to pick the distribution
of points on the grid. For both types of problems, it will also always
be used in calls to edge generators, separation routines,
etc. Negative arguments, or an argument of zero, will result in the
current time being used.

By default, Camargue will do a loop of pivoting and cutting for as
long as possible, and then begin a so-called Augment-Branch-Cut (ABC)
search. The flag option `-P` will disable branching, attempting a "pure"
primal cutting plane solution method instead.

Camargue implements several different node selection rules for guiding
the ABC search -- these can be specified by passing options to
`-b`. For example, `-b 3` will run a depth-first search traversal of
the ABC tree, whereas `-b 2` will do a primal variant of the familiar
best-bound, or best-first, search. The default option is specified
with `-b 0`: a best-bound search will be interleaved into a so-called
best-tour search. A pure best-tour search is selected with `-b 1`.

The flag option `-S` is available to specify sparse solution
mode. In this mode, Camargue will run no edge pricing of any kind; it
will just generate an initial edge set as above and try to prove that
a given starting tour is optimal for this edge set. This option is
required for the use of primal safe Gomory cut separation, described
in the [external dependencies](externals/extdeps.md).

Moreover, users can also specify cut generation style with
`-c`. Camargue contains implementations of certain primal separation
algorithms based on the research of Letchford and Lodi, as well as
Fleischer, Letchford and Lodi. The option `-c 0` will select these
algorithms, as well as heuristics for fast blossoms and block
combs. The `-c 1` option, the default, will supplement these with some
more exotic standard separation routines implemented in the Concorde
TSP solver, such as cut tightening, double deckers, path inequalities,
comb teething, and local cuts. Camargue will struggle to solve most
instances with the `-c 0` option, but it can still be used.

Finally, there is the `-l [int or float]` option to specify a target
lower bound for the solver. Since Camargue works by trying to augment
starting tours, you may be satisfied with terminating the solver
prematurely if some target objective value is met. Values supplied
will be rounded up to their integer ceiling, with the interpretation
that integer values represent tour lengths and floating points
represent dual lower bounds.
