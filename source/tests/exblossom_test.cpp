#include "config.hpp"

#ifdef CMR_DO_TESTS


#include "datagroups.hpp"
#include "process_cuts.hpp"
#include "blossoms.hpp"
#include "separator.hpp"
#include "err_util.hpp"

#include <catch.hpp>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

SCENARIO ("Tiny primal exact blossom separation",
          "[Sep][ExBlossoms][tiny]") {
    using namespace CMR;

    GIVEN ("The tiny instance blossom6 with an obvious blossom") {
        THEN ("Blossoms are found") {
            Graph::CoreGraph core_graph;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            vector<double> lp_edges;

            Data::make_cut_test("test_data/blossom6.tsp",
                                "test_data/tours/blossom6.sol",
                                "test_data/subtour_lp/blossom6.sub.x",
                                core_graph, b_dat, lp_edges, s_dat);
            Sep::CutQueue<Sep::ex_blossom> blossom_q;
            LP::ActiveTour act_tour(core_graph, b_dat);
            Sep::ExBlossoms ex_b(core_graph.get_edges(),
                                 act_tour, s_dat,
                                 blossom_q);
            REQUIRE(ex_b.find_cuts());
        }
    }
}

SCENARIO ("Primal exact blossom separation",
          "[Sep][ExBlossoms]") {
    using namespace CMR;
    vector<string> probs {
        "pr76",
        // "d493",
        // "att532",
        // "pr1002",
        // "rl1304",
        // "d2103",
        // "pcb3038",
    };

    for (string &prob : probs) {
        GIVEN ("LP vectors and tours for " + prob) {
            string
            probfile = "problems/" + prob + ".tsp",
            solfile = "test_data/tours/" + prob + ".sol",
            blossomfile = "test_data/blossom_lp/" + prob + ".2m.x",
            subtourfile = "test_data/subtour_lp/" + prob + ".sub.x";

            Graph::CoreGraph core_graph;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            vector<double> lp_edges;

            WHEN ("The solution is in the blossom polytope") {
                Data::make_cut_test(probfile, solfile, blossomfile,
                                    core_graph, b_dat, lp_edges, s_dat);

                THEN ("No blossoms are found") {
                    Sep::CutQueue<Sep::ex_blossom> blossom_q;
                    LP::ActiveTour act_tour(core_graph, b_dat);
                    Sep::ExBlossoms ex_b(core_graph.get_edges(),
                                         act_tour, s_dat,
                                         blossom_q);

                    REQUIRE_FALSE(ex_b.find_cuts());
                }
            }

            WHEN ("The solution is in the subtour polytope") {
                Data::make_cut_test(probfile, solfile, subtourfile,
                                    core_graph, b_dat, lp_edges, s_dat);

                THEN ("ExBlossoms succeeds when FastBlossoms does") {
                    Sep::CutQueue<Sep::ex_blossom> blossom_q;
                    LP::ActiveTour act_tour(core_graph, b_dat);
                    Sep::ExBlossoms ex_b(core_graph.get_edges(),
                                         act_tour, s_dat,
                                         blossom_q);

                    int found_ex = ex_b.find_cuts();

                    Data::KarpPartition kpart;

                    Sep::Separator sep(core_graph.get_edges(), act_tour, s_dat,
                                       kpart, 99);

                    int found_fast = sep.fast2m_sep();
                    INFO("Known issue: ExBlossoms misses some fast blossoms");
                    INFO("Probably need to implement a work_blossom function");
                    CHECK(found_ex >= found_fast);

                    if (found_ex < found_fast) {
                        auto &f2m_q =
                        const_cast<Sep::LPcutList &>(sep.fastblossom_q());
                        cout << f2m_q.size() << " fast2m cuts\n";
                        for (auto it = f2m_q.begin(); it; it = it->next) {
                            CCtsp_print_lpcut_in(it);
                        }
                    }
                }
            }
        }
    }
}

SCENARIO ("Black box ExBlossoms testing",
          "[.Sep][.ExBlossoms][cut_ecap]") {
    using namespace CMR;
    vector<string> probs {
        "pr76",
        "d493",
        "att532",
        "pr1002",
        "rl1304",
        "d2103",
        "pcb3038",
    };

    for (string &prob : probs) {
        GIVEN ("LP vectors and tours for " + prob) {
            string
            probfile = "problems/" + prob + ".tsp",
            solfile = "test_data/tours/" + prob + ".sol",
            blossomfile = "test_data/blossom_lp/" + prob + ".2m.x",
            subtourfile = "test_data/subtour_lp/" + prob + ".sub.x";

            Graph::CoreGraph core_graph;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            vector<double> lp_edges;

            Data::make_cut_test(probfile, solfile, subtourfile,
                                core_graph, b_dat, lp_edges, s_dat);

            WHEN ("We set up ecaps as in find_cuts") {
                vector<int> &sup_inds = s_dat.support_indices;
                vector<double> &sup_ecap = s_dat.support_ecap;
                vector<int> &sup_elist = s_dat.support_elist;
                const vector<int> &tour_edges = b_dat.best_tour_edges;

                vector<double> cut_ecap;

                cut_ecap = sup_ecap;

                for (auto i = 0; i < sup_inds.size(); ++i) {
                    int index = sup_inds[i];
                    if (tour_edges[index] == 1)
                        cut_ecap[i] = 1 - sup_ecap[i];
                }

                vector<double> saved_copy = cut_ecap;

                THEN ("cut_ecap is restored after emulating the main loop") {
                    for (auto i = 0; i < sup_inds.size(); ++i) {
                        int cut_ind = sup_inds[i];
                        int tour_entry = tour_edges[cut_ind];
                        double orig_weight = 1.0;
                        double changed_weight = 1.0;

                        if (tour_entry == 0) {
                            orig_weight = sup_ecap[i];
                            changed_weight = 1 - sup_ecap[i];
                        } else if (tour_entry == 1) {
                            orig_weight = 1 - sup_ecap[i];
                            changed_weight = sup_ecap[i];
                        }

                        cut_ecap[i] = changed_weight;

                        //reverts cut_ecap when it goes out of scope.
                        auto ecap_guard = util::make_guard([&cut_ecap, i,
                                                            orig_weight]
                                                           { cut_ecap[i] =
                                                               orig_weight; }
                            );
                    }

                    REQUIRE(cut_ecap == saved_copy);
                }

            }

        }
    }
}
#endif
