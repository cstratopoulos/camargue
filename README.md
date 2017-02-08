Camargue	{#mainpage}
========

This is the README for Camargue, a TSP solver based on primal
cutting plane methods. Camargue tries to move from tour to tour, either
proving that a tour is optimal or finding a better one. This page contains
information on installing and using Camargue. 

See below
for installation and usage. To get a sense of the layout of the code, a good
starting point would be the documentation for CMR::Solver.

Installation
------------

For users building the code on one of the UW linux.math environments,
from the camargue main directory just do

    cp scripts/Makefile.template.UWlinux Makefile

The compiler and CPLEX directories are already specified, so you just
need to do the steps below concerning symlinks to Concorde.

Camargue is a C++11 application, so most importantly you will need a
compiler that is compliant with most of the C++11 standard. In fact the code
has been tested with GCC 4.7, which technically does *not* support the
C++11 standard in full. Anything after 4.7 or later should be fine,
and `g++` on a Mac (which is actually an alias of Apple `clang`)
should be fine too. Your compiler should be specified by editing the
`CC` definition in the `Makefile`. A bit more detail and example
options are given there. 

Camargue relies heavily on two main external dependencies:
- The TSP solver/combinatorial optimization library
[Concorde](http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm).

- The linear programming solver [CPLEX
12](http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud)

Both must be installed before proceeding.

After installing Concorde, go into the directory `camargue/externals`
and create a symlink to the `concorde` directory. That is, `concorde`
must be the name of the folder containing the files `TSP`, `INCLUDE`,
`CUT`, etc. 

After installing CPLEX, open the `Makefile` and edit the definitions
`CPXDIR` and `CPX_LIB`. Details and examples are given in the
`Makefile`.

No further edits to the Makefile should be necessary. After that, run
`camargue/configure.sh`. The most important task done by this script
is to edit the file `concorde/INCLUDE/tsp.h` to remove instances of
the keyword `new` appearing in function prototypes.

The script also attempts to configure Camargue to use certain
external dependencies. These are all discussed in detail in @ref
extdeps. With none of these present, or if the configuration script
fails at this step somehow, the project should still compile. To see
if everything went OK, you can check the file config.hpp. If you went
with a barebones install and later added external enhancements, you
should be able to use them with `make clean` and then by running
`configure.sh` again. 

After performing all these steps, you should be able to compile and
run Camargue by running `make` from the main directory, creating the
`camargue` executable.

Additionally, if you like, you can run the unit
tests, benchmarks, and experiments that I used to develop
Camargue (and write my thesis!) by using the recipe `make test`. This
requires that the Catch unit testing framework be
installed. Information on this is given in @ref extdeps, and specific
information on invoking the unit tests is in @ref unittests.

Usage
------

This heading is about standard command line usage of Camargue. For
information on running tests/benchmarks, see @ref unittests.

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
of solution files supported is *not* the TSPLIB `.tour` format. Rather,
it should be a file whose first line is the instance node count, with
the following lines (with arbitrary spacing/indentation) giving
zero-indexed ordering of the nodes. Camargue will make sure no obvious
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
code.

Also for both styles of problems, you can pass a random seed with
`-s`. This is to allow reproducibility through all areas of the
code. For a random problem, this will be used to pick the distribution
of points on the grid. For both types of problems, it will also always
be used in calls to edge generators, separation routines,
etc. Negative arguments, or an argument of zero, will result in the
current time being used.

Finally, the flag option `-S` is available to specify sparse solution
mode. In this mode, Camargue will run no edge pricing of any kind; it
will just generate an initial edge set as above and try to prove that
a given starting tour is optimal for this edge set. This option is
required for the use of primal safe Gomory cut separation; see @ref
extdeps for more info.

@TODO Should be posisble to specify starting edge sets from file for
TSPLIB instances. 