#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "core_lp.hpp"
#include "datagroups.hpp"
#include "separator.hpp"
#include "util.hpp"

#include <array>
#include <vector>
#include <string>
#include <iostream>

#include <catch.hpp>

using std::array;
using std::vector;

using std::string;
using std::cout;

SCENARIO ("Pivoting and adding cuts",
         "[LP][Sep][CoreLP][Separator]") {
    using namespace CMR;
    vector<string> probs{"dantzig42", "st70", "pr76", "lin105",
                         "lin318", "d493", "att532", "pr1002", "rl1304",
                         "d2103", "pr2392"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("We pivot to a new solution") {
                THEN ("We can find cuts, add them, and pivot back") {
                    Data::Instance inst("problems/" + prob + ".tsp", 99);
                    Graph::CoreGraph core_graph(inst);
                    Data::BestGroup b_dat(inst, core_graph);
                    LP::CoreLP core(core_graph, b_dat);
                    Data::KarpPartition kpart(inst);

                    vector<double> tour1 = core.lp_vec();

                    LP::PivType piv = core.primal_pivot();
                    cout << "Pivoted type: " << piv
                         << ", piv val: " << core.get_objval() << "\n";

                    vector<double> piv1 = core.lp_vec();

                    vector<int> island;
                    Data::SupportGroup s_dat(core_graph.get_edges(),
                                             piv1, island,
                                             inst.node_count());

                    Sep::Separator control(core_graph.get_edges(),
                                           core.get_active_tour(), s_dat,
                                           kpart, 99);

                    bool fast2m = control.fast2m_sep();
                    bool blkcomb = control.blkcomb_sep();
                    bool dp = control.simpleDP_sep();
                    bool con = control.connect_sep();
                    bool found = fast2m || blkcomb || dp || con;
                    REQUIRE(found);

                    REQUIRE_NOTHROW(core.pivot_back(false));

                    vector<double> tour2 = core.lp_vec();

                    REQUIRE(tour1 == tour2);


                    std::array<Sep::LPcutList*,
                               4> qlist{
                        &control.segment_q(),
                        &control.fastblossom_q(),
                        &control.blockcomb_q(),
                        &control.connect_cuts_q()};

                    for (auto qptr : qlist)
                        REQUIRE_NOTHROW(core.add_cuts(*qptr));

                    REQUIRE_NOTHROW(core.add_cuts(control.simpleDP_q()));

                    piv = core.primal_pivot();
                    cout << "Pivoted type: " << piv
                         << ", piv val: " << core.get_objval() << "\n";

                    vector<double> piv2 = core.lp_vec();

                    REQUIRE(piv1 != piv2);


                    cout << "\n\n";
                }
            }
        }
    }
}

#endif //CMR_DO_TESTS
