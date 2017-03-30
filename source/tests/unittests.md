Catch Unit Tests
=================

This file describes how to interactively use the Catch unit tests that were
used to develop Camargue. For information on installing Catch and
building in test mode, see the [general external dependencies
documentation](../../externals/extdeps.md). Henceforth we'll assume
Catch is installed, and that all the test data has been generated with
a call to

    ./make_test_data.sh

Running the tests
------------------

Compiling with `make test` defines a new main function for Camargue
that is provided by Catch. You can run

    ./camargue

but this might not be what you want: this runs *every single* test
case programmed, and will probably take 20-40 minutes to run. More
likely, you want to examine tests for smaller chunks of code.

Catch offers a tagging system to label test cases, and I have used
this to indicate logical units of the code that
are covered by a given test. These mostly correspond to namespaces,
class names, or function names. To give some examples,

    ./camargue '[Sep]'

will run all tests associated to the namespace CMR::Sep, for matters
related to separation routines. And

    ./camargue '[HyperGraph]'

will run tests corresponding to the class CMR::Sep::HyperGraph. To
increase the specificity even further,

    ./camargue '[get_coeffs]'

will test the method CMR::Sep::HyperGraph::get_coeffs. In some special
cases, the tag will refer to a class member. For example, I used the
tests in

    ./camargue '[adj_zones]'

to do brute force testing of the adjacency zone technique I developed
for candidate teeth separation.

Catch also offers a natural syntax for combining or excluding test
cases. The logical "and" syntax looks like this:

    ./camargue '[Sep][LP]'

will run all tests tagged with both the namespaces CMR::Sep and
CMR::LP, indicating that they should test operations that make
nontrivial use of facilities from both namespaces. Comma separated
lists such as `[Sep],[LP]` indicate a logical "or", and a tilde can
exclude a tag. For example

    ./camargue '[Solver]~[random]'

will run tests of CMR::Solver while excluding the tag `[random]`,
which I used to label test cases with randomly generated Euclidean
instances. You can try typing in any names of namespaces and classes
that you like to see if tests exist for them, but you can also run

    ./camargue -t

to list all the test case tags.

In most of these cases, Camargue and Catch will probably report that
no errors occurred, but it might not be clear what was being done in
the first place. By default, Catch only reports names and conditions
of test cases when they fail. To override this behavior, run the tests
with the flag `-s`. Then, every test case will explicitly report what
was being tested along with assertions that succeeded or failed. So

    ./camargue '[DPwitness]' -s

will run tests of the class CMR::Sep::DPwitness, and each case will be
reported like

> Scenario: Finding simple DP inequalities via karp partition witnesses
>
>    Given: A karp partition and candidate teeth for pcb3038
>
>     Then: We can get simple DP inequalities in a mini cutgraph


This is probably more than enough info to get started running the
Camargue test cases, but if you want to know more about using tags, or
how Catch works, the documentation on the [Catch GitHub
page](https://github.com/philsquared/Catch) is helpful.

Special tests and tags
-----------------------

In truth, running `./camargue` won't run quite every single test:
Catch also offers the tag `[!hide]` or `[.]` to mark test cases as hidden, so
that they don't run by default. I have used the hiding tag in only a
couple situations.

Tests with the tag `[valgrind]` are marked as hidden, but can be
invoked directly with

    ./camargue '[valgrind]'

These tests mostly consist of cases where I have deliberately caused
constructor errors in classes which manage C resources and are
responsible for freeing managed memory upon destruction. As the name
suggests, you are invited to run them with
[Valgrind](http://valgrind.org/) to see that no memory is leaked
during exception conditions.

Similarly, tests with the tag `[!shouldfail]` are marked by Catch as
hidden by default. As the name suggests, these are tests of the error
handling control flow in Camargue. They ensure that the program
terminates as expected if provided with bad inputs.

Another notable use of the hidden tags is described in the next section.

Thesis Experiments
--------------------

In developing Camargue, I have used Catch not only to test the code
for bugs, but to generate benchmarks, tables, and figures that are
reported in my thesis. These are implemented through the hidden tag
system, as described in the previous section. In particular you can
try tests like

    ./camargue '[.benchmark]'
    ./camargue '[.figure]'
    ./camargue '[.table]'

to run tests that correspond to thesis results. Tests with
`'[figure]'` will usually write tours/coordinates/LP solutions to
file, which can then be rendered with some sort of graph visualization
software. With `'[.benchmark]'` or `'[.table]'`, sometimes the tests
will produce almost verbatim a tab-separated table which was
then copied into my thesis. In other cases, a bit of extra processing
was used to format or interpret the data. For all of these, I would
suggest running with the `-s` flag, as described previously, to get an
idea of what the results actually mean. Additionally, you may want to
browse the source of the test case as well to determine which summary
statistic was being used.

Note the `.` prefix on all the test cases: they are all hidden by
default because many of them run *extremely* slowly, given that they
try to provide empirical evidence for one approach being faster than
another. On top of this, some of them run these slow trials 5 or 10 or
20+ times, so as to record mean CPU times.

At the end of the chapter on simple DP cuts, I report on several
different experiments with the TSPLIB instance pla85900. A bit more
work is required with these, since they involve modifying parameters
which are determined at compile time. They are associated to the tags
`'[.sdp-pla85]'` and such, and more detailed instructions for running
them are given in `simpleDP_test.cpp`.

The tag `'[experiment]'` is also employed in a very limited fashion:
tests with this tag were used to get a sense of the workings of
certain aspects of the solution process, so as to guide implementation
choices.

Reading the tests
-------------------------

Hopefully the use of tags and the `-s` flag described above will give
a satisfying impression of what the tests are doing while making it
easy to target specific parts of the code for testing. If some of
these are unclear (or if you're unconvinced about test successes!),
you are invited to examine the test implementation files
yourself. Test cases are written in the Behavior-Driven Development
style, usually with a full sentence describing the behavior of each
part. The hope is that this, combined with the natural syntax of
Catch, will make the tests more or less self-documenting.

That being said, the source files themselves are much less readable
than the rest of the source code. There are around 4,500 source lines of test
code, and a lot of it is repetitive boilerplate that could probably be
abstracted somehow. Similarly, some of the tests that produce
formatted tables do so using some fairly lamentable array-indexing
tricks.