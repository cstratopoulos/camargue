#include "config.hpp"

#ifdef CMR_DO_TESTS


#include "datagroups.hpp"
#include "process_cuts.hpp"
#include "blossoms.hpp"
#include "separator.hpp"

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
            Data::GraphGroup g_dat;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            vector<double> lp_edges;

            Data::make_cut_test("test_data/blossom6.tsp",
                                "test_data/tours/blossom6.sol",
                                "test_data/subtour_lp/blossom6.sub.x",
                                g_dat, b_dat, lp_edges, s_dat);
            Sep::CutQueue<Sep::ex_blossom> blossom_q;
                    
            Sep::ExBlossoms ex_b(g_dat.core_graph.get_edges(),
                                 lp_edges, b_dat, s_dat,
                                 blossom_q);
            REQUIRE(ex_b.find_cuts());
        }
    }
}

SCENARIO ("Primal exact blossom separation",
          "[Sep][ExBlossoms]") {
    using namespace CMR;    
    vector<string> probs {
        "att532",
        "d2103",
        "d493",
        "pcb3038",
        "pr1002",
        "pr76",
        "rl1304",
    };

    for (string &prob : probs) {
        GIVEN ("LP vectors and tours for " + prob) {
            string
            probfile = "problems/" + prob + ".tsp",
            solfile = "test_data/tours/" + prob + ".sol",
            blossomfile = "test_data/blossom_lp/" + prob + ".2m.x",
            subtourfile = "test_data/subtour_lp/" + prob + ".sub.x";

            Data::GraphGroup g_dat;
            Data::BestGroup b_dat;
            Data::SupportGroup s_dat;
            vector<double> lp_edges;
            
            WHEN ("The solution is in the blossom polytope") {
                Data::make_cut_test(probfile, solfile, blossomfile,
                                    g_dat, b_dat, lp_edges, s_dat);

                THEN ("No blossoms are found") {
                    Sep::CutQueue<Sep::ex_blossom> blossom_q;
                    
                    Sep::ExBlossoms ex_b(g_dat.core_graph.get_edges(),
                                         lp_edges, b_dat, s_dat,
                                         blossom_q);

                    REQUIRE_FALSE(ex_b.find_cuts());
                }
            }

            WHEN ("The solution is in the subtour polytope") {
                Data::make_cut_test(probfile, solfile, subtourfile,
                                    g_dat, b_dat, lp_edges, s_dat);

                THEN ("ExBlossoms succeeds when FastBlossoms does") {
                    Sep::CutQueue<Sep::ex_blossom> blossom_q;
                    
                    Sep::ExBlossoms ex_b(g_dat.core_graph.get_edges(),
                                         lp_edges, b_dat, s_dat,
                                         blossom_q);

                    int found_ex = ex_b.find_cuts();

                    Graph::TourGraph TG(b_dat.best_tour_edges,
                                        g_dat.core_graph.get_edges(),
                                        b_dat.perm);
                    Data::KarpPartition kpart;

                    Sep::Separator sep(g_dat, b_dat, s_dat, kpart, TG);

                    int found_fast = sep.fast2m_sep();

                    REQUIRE(found_ex >= found_fast);
                }
            }
        }
    }
}
#endif
