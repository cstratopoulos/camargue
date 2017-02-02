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
[Catch](https://github.com/philsquared/Catch). Catch is a header-only,
C++-based application. To use it with Camargue, just put the [Catch
header
file](https://raw.githubusercontent.com/philsquared/Catch/master/include/catch.hpp)
(or a symlink to it) in the 'camargue/externals' folder.

To compile Camargue in testing mode, use the recipe `make test`. This
edits a line in config.hpp to enable testing, `make`s the project,
then reverts config.hpp afterward. Thus, subsequent calls to `make`
should just build the normal command line executable. This will *not*
happen if, for some reason, the compilation process is killed midway
through `make test`. In this case, `make remake` should revert
config.hpp and allow the normal Camargue main to be built. 


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