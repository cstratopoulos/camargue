External Dependencies	{#extdeps}
=====================

Camargue  would not exist without heavy reuse of facilities
from Concorde, CPLEX, and C++11/the C++ STL. The philosophy has been
to reuse code from these sources wherever reasonably (and sometimes
unreasonably) possible. In light of this, remaining dependencies have
been implemented on an as-you-like basis.

In all cases, it should be sufficient to `make clean` and then run
`configure.sh` again if a new dependency has been added. 

Safe Gomory Cuts
----------------

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

@TODO The raw download will not compile. Maybe provide a script which automates
the necessary changes, or provide easy access to the modified version.

Create a softlink to the directory
`safemir` in the `camargue/externals` folder. So `safemir` must be the
name of the file containing the folder `src`.

Catch Unit Tests
----------------

This section just describes installation of the Catch unit testing
framework. For info on running the tests, and suggested usage, see
@ref unittests. 

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
should just build the normal command line executable. This will *not*
happen if, for some reason, the compilation process is killed midway
through `make test`. In this case, `make remake` should revert
config.hpp and allow the normal Camargue main to be built.

If you want to run the tests, there are some requirements on
what your `camargue` directory has to look like. Most importantly, all
the tests search for TSPLIB instances in a folder called `problems` in
the Camargue directory. The 
largest problem used in any of the tests is `usa13509`. Thus, if you
have a folder containing the TSPLIB instances from
[here](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/),
and you have extracted all problems of size at most as big as
`usa13509`, just create a symlink to the containing folder in the
`camargue/` directory. The other requirement is the `test_data`
folder. This folder contains three subdirectories:

- `tours` for (often suboptimal, sometimes deliberately terrible) TSP tours
- `subtour_lp` for files with names like `att532.sub.x`
- `blossom_lp` for files with names like `d493.2m.x`

The latter two are cases where Concorde was used with the flag `-i` or
`-I` to solve over the blossom or subtour polytope, respectively. If
you don't have these files, the script TODO can be used to generate
them, provided

- you have linked to a working build of Concorde in `externals`, and
- you have created the `problems` symlink in the camargue directory.

The exception to these are a family of extremely tiny examples that I
created, `blossom6.tsp`, `comb9.tsp`, `fleisA9.tsp` and
`fleisB9.tsp`. The instances `blossom6` and `comb9` are familiar
blackboard examples used to illustrate subtour and comb
inequalities. The instances `fleisA9` and `fleisB9` are an attempt to
recreate Figures 2 and 5 from Fleischer et al.'s 2006 paper on simple
DP inequalities. These TSP instances have been placed in the
`test_data` directory, with tours and LP solutions in subdirectories
as appropriate. 

If you want to develop or modify tests, repeatedly calling `make test`
will cause slower compile times, so it would be better to manually
change the `CMR_DO_TESTS` line in config.hpp and change it back when
you are done. 


OpenMP Parallelism
-------------------

There are a few separation routines in the code that can be
implemented to run in parallel. I have chosen to implement this using
the [OpenMP](http://www.openmp.org/) standard for an extremely simple
approach. Unlike the other external dependencies, this one is not a
download. Rather, the `configure.sh` script will try to test if your
compiler supports OpenMP (OMP), and edit the `Makefile`
accordingly. 

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

Sorting with Timsort
---------------------

Timsort is an enhancement applied during candidate tooth
generation. In my thesis I show that the subroutine I use finds lists
of teeth that are often already sorted, or consist of a few sorted
subsequences. Compared with quicksort, Timsort is an algorithm that
performs much better on already-sorted or partially-sorted data. But
again, the sorting accounts for an extremely small portion of the
computation in finding candidate teeth.

I have chosen to invoke Timsort using the implementation by Fuji Goro,
available [here](https://github.com/gfx/cpp-TimSort). To use it, just
put the
[header](https://raw.githubusercontent.com/gfx/cpp-TimSort/master/timsort.hpp)
in the `externals` folder. If `configure.sh` detects that it is not
present, the code will just use the C++ STL `std::sort`.