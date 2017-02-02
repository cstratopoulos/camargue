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

To use it, follow the instructions in the `Makefile` to create softlinks
to the CPLEX root directory and to the folder containing
`libcplex.a`. After that, run `make` for good measure just to ensure
everything is working. Camargue itself will compile only the
headers/objects it needs, using a bit of a crude approach where the
needed header files (and certain implementation files) are included
directly. Thus no edits to the Camargue `Makefile` should be
necessary.

If the installation looks OK, create a softlink to the directory
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


