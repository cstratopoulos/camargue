External Dependencies
=====================

Camargue  would not exist without heavy reuse of facilities
from Concorde, CPLEX, and C++11/the C++ STL. The philosophy has been
to reuse code from these sources wherever reasonably (and sometimes
unreasonably) possible. In light of this, remaining dependencies have
been implemented on an as-you-like basis.

As mentioned in the README,

    ./cmr_install.py --full

does a full install with all the dependencies described below. To add
individual ones, it should be sufficient to run `make clean` and then
pass additional flags to the install script with the `-W` flag. The
add-ons, and their corresponding flags, are as follows:

- Primal separation of safe Gomory cuts, `-W safemir`
- Unit tests and benchmarks with Catch, `-W catch`
- Shared-memory parallelism with OpenMP, `-W omp`

Safe Gomory Cuts
-------------------------

In "Primal Cutting Plane Algorithms Revisted", Letchford and Lodi
describe how to generate Gomory cuts that solve the primal separation
problem. To avoid issues with edge pricing, Gomory cuts are only used
on sparse instances. In these cases, though, their effect on the
solution process is quite dramatic. For an idea of how they are called
in Camargue, see the documentation for CMR::Sep::SafeGomory.

Safe Gomory cuts are implemented
using the code of William Cook, Sanjeeb Dash, Ricardo Fukasawa, and
Marcos Goycoolea, provided as a supplement to their paper "Numerically
Safe Gomory Mixed-Integer Cuts" in the Informs Journal of
Computing. The code is available as a supplement from the publisher
[here](https://www.informs.org/Pubs/IJOC/Online-Supplements/Volume-21-2009/Cook-Dash-Fukasawa-Goycoolea).

If requested, the Camargue install script will downloaded the gzipped
tarball at the link above, extracting it to `camargue/externals`. The
script will then make some minor edits to some files in
`safemir/src`. Some of the edits just prevent compiler warnings, but
others are required for the project to compile.

If you have already downloaded the code with the link above, create a
symlink to the directory `safemir` in the `camargue/externals`
folder. So `safemir` must be the name of the file containing the
folder `src`. If you already had the code downloaded, you will still
need to make the automated changes just described. The Camargue
install script will check if the folder is already present, in which
case it will only make the edits, rather than re-downloading and extracting.

Catch Unit Tests
---------------------------

This section just describes installation of the Catch unit testing
framework. For info on running the tests, and suggested usage, see
[here](../source/tests/unittests.md).

The current version of Camargue has been developed in a test driven
fashion with the help of the unit testing framework
[Catch](https://github.com/philsquared/Catch): C++ Automated Test
Cases in Headers. To use it with Camargue, just put the
[Catch
header](https://raw.githubusercontent.com/philsquared/Catch/master/include/catch.hpp)
(or a symlink to it) in the 'camargue/externals' folder.

To compile Camargue in testing mode, use the recipe `make test`. This
edits a line in config.hpp to enable testing, `make`s the project,
then reverts config.hpp afterward. Thus, subsequent calls to `make`
should just build the normal command line executable. This will _not_
happen if, for some reason, the compilation process is killed midway
through `make test`. In this case, `make remake` or `make undeftest
all` will revert `config.hpp`, again allowing the Camargue solver to
be built.

If you want to run the tests, there are some requirements on
what your `camargue` directory has to look like. Most importantly, all
the tests search for TSPLIB instances in a folder called `problems` in
the Camargue directory, so at `camargue/problems`. Thus, if you
have a folder containing the TSPLIB instances from
[here](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/),
and you have extracted the problems, just create a symlink to the
containing folder in the `camargue/` directory. The other requirement
is the `test_data`
folder. This folder contains three subdirectories:

- `tours` for (often suboptimal, sometimes deliberately terrible) TSP tours
- `subtour_lp` for files with names like `att532.sub.x`
- `blossom_lp` for files with names like `d493.2m.x`

The latter two are cases where Concorde was used with the flag `-i` or
`-I` to solve over the blossom or subtour polytope, respectively.
Running `scripts/make_test_data.py` will generate these files for you.

- you have linked to a working build of Concorde in `externals`, and
- you have created the `problems` symlink in the camargue directory.

If you have some files but not others, the script will do a simple
check before actually invoking Concorde.

The tests used to develop Camargue, and to write my thesis, contain
examples of size up to and including the largest TSPLIB instance,
`pla85900`. Although Concorde can very quickly compute subtour
polytope optima for such instances, you may want to skip them for
quicker generation of the majority of test cases. By default,
`make_test_data.py` skips examples of size bigger than 15,000 nodes,
but you can override this behaviour with the flag `--large`.

The exception to the cases described are a family of extremely tiny
examples that I created, `blossom6.tsp`, `comb9.tsp`, `fleisA9.tsp` and
`fleisB9.tsp`. The instances `blossom6` and `comb9` are familiar
blackboard examples used to illustrate subtour and comb
inequalities. The instances `fleisA9` and `fleisB9` are an attempt to
recreate Figures 2 and 5 from Fleischer et al.'s 2006 paper on simple
DP inequalities. These TSP instances have been placed in the
`test_data` directory, with tours and LP solutions in subdirectories
as appropriate.

If you want to develop or modify tests (or do the special tests with
`pla85900`), `make test` will cause slow compile times, and it will
not be possible to manually enable OpenMP, if desired. For this, I
have provided the recipe `make develop_test`, which modifies
`config.hpp` to build tests, without changing it back after. To go
back to normal from this, you can use `make undeftest all` or `make remake`.


OpenMP Parallelism
---------------------------

There are a few separation routines in the code that can be
implemented to run in parallel. I have chosen to implement this using
the [OpenMP](http://www.openmp.org/) standard for an extremely simple
approach. Unlike the other external dependencies, this one is not a
download. Rather, if requested, `cmr_install.py` will try to test if
your compiler supports OpenMP (OMP), editing the `Makefile`
accordingly. Unlike the other examples, there is no script that will
download anything: if your compiler doesn't support OMP it is probably
not worth it to download and install a new one just for that.

In "Primal Separation Algorithms", Letchford and Lodi observe that the
minimum cut computations in their blossom separation algorithm are
"independent of each other", so I have implemented these to run
concurrently. And in my approach to primal simple DP separation, it is
possible to search the Karp partitioned witness cutgraphs in
parallel. From a practical point of view in either case the speedup is
pleasant but not earth shattering, as neither separation routine is
called terribly often.

The OMP standard dictates that if the compiler does not support
OMP `#pragma`s, they are simply ignored and the result is still valid
code. However in my implementations there is a bit of added overhead
for memory management or error checking, so the result is not as clean
as the implementation that could be used in the serial case.
