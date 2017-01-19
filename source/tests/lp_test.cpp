#include "config.hpp"

#ifdef CMR_DO_TESTS

#include <cplex.h>

#include "lp_interface.hpp"
#include "core_lp.hpp"
#include "separator.hpp"

#include "err_util.hpp"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <catch.hpp>

using std::cout;
using std::vector;
using std::unique_ptr;
using std::string;
using std::pair;

SCENARIO ("Performing single pivots",
          "[LP][CoreLP][primal_pivot][Sep][Separator]") {
    vector<string> probs{"pr76", "a280", "p654", "pr1002", "rl1304"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("We pivot to a new solution") {
                THEN ("The LP vector changes") {
                    CMR::Data::Instance inst("problems/" + prob + ".tsp", 99);
                    CMR::Data::GraphGroup g_dat(inst);
                    CMR::Data::BestGroup b_dat(inst, g_dat);
                    CMR::LP::CoreLP core(g_dat, b_dat);

                    vector<double> tourx = core.lp_vec();
                    double tourlen = core.get_objval();

                    REQUIRE_NOTHROW(core.primal_pivot());

                    vector<double> pivx = core.lp_vec();
                    double pval = core.get_objval();

                    REQUIRE(pval != tourlen);
                    REQUIRE(tourx != pivx);
                    cout << "\tPivot val: " << pval << "\n";

                    AND_WHEN("We pivot back") {
                        THEN("The tour vector changes back.") {
                            REQUIRE_NOTHROW(core.pivot_back());
                            REQUIRE(core.get_objval() == tourlen);
                            REQUIRE(tourx == core.lp_vec());

                            CMR::Data::SupportGroup s_dat;

                            s_dat.reset(inst.node_count(),
                                        g_dat.core_graph.get_edges(), pivx,
                                        g_dat.island);

                            CMR::Graph::TourGraph TG(b_dat.best_tour_edges,
                                                     g_dat.core_graph.
                                                     get_edges(), b_dat.perm);
                            CMR::Data::KarpPartition kpart;
                            CMR::Sep::Separator sep(g_dat, b_dat, s_dat,
                                                    kpart);

                            bool found_some = false;

                            if (sep.segment_sep(TG)) {
                                found_some = true;
                                core.add_cuts(sep.segment_q());
                            } else if (sep.fast2m_sep(TG)) {
                                found_some = true;
                                core.add_cuts(sep.fastblossom_q());
                            } else if (sep.blkcomb_sep(TG)) {
                                found_some = true;
                                core.add_cuts(sep.blockcomb_q());
                            } else if (sep.connect_sep(TG)) {
                                found_some = true;
                                core.add_cuts(sep.connect_cuts_q());
                            }

                            REQUIRE(found_some);
                            
                            AND_WHEN("We pivot again after adding cuts") {
                                THEN("The tour vector changes again.") {
                                    core.primal_pivot();
                                    vector<double> piv2x = core.lp_vec();
                                    double piv2val = core.get_objval();
                                    cout << "\tPiv2 val: " << piv2val << "\n";

                                    REQUIRE(piv2x != pivx);
                                    REQUIRE(piv2val != pval);
                                    cout << "\t Pivot delta: "
                                         << (piv2val - pval) << "\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

SCENARIO ("Consructing a Core LP",
          "[LP][CoreLP]") {
    vector<string> probs{"dantzig42", "lin318", "d493", "pr1002", "pcb3038"};

    for (string& prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("A core LP is constructed") {
                CMR::Data::Instance inst("problems/" + prob + ".tsp", 99);
                CMR::Data::GraphGroup g_dat(inst);
                CMR::Data::BestGroup b_dat(inst, g_dat);
                THEN ("Its constructor doesn't throw.") {
                    REQUIRE_NOTHROW(CMR::LP::CoreLP core(g_dat, b_dat));
                }

                AND_THEN ("The degree LP is feasible at the best tour") {
                    CMR::LP::CoreLP core(g_dat, b_dat);
                    vector<double> feas;
                    vector<double> tour = core.lp_vec();
                    REQUIRE_NOTHROW(core.get_row_infeas(tour, feas, 0,
                                                        core.num_rows() - 1));
                    bool found_infeas = false;
                    for(double &stat : feas) {
                        if (stat) {
                            found_infeas = true;
                            break;
                        }
                    }

                    REQUIRE_FALSE(found_infeas);
                }
            }
        }
    }
}

SCENARIO ("Constructing LP Relaxations",
          "[.LP][.Relaxation][valgrind]") {
    WHEN("A relaxation is constructed"){
        THEN("Its constructor doesn't throw"){
            REQUIRE_NOTHROW(CMR::LP::Relaxation rel);            
        }
    }

    GIVEN ("A constructed Relaxation") {
        CMR::LP::Relaxation rel;
        WHEN ("Another one is constructed") {
            CMR::LP::Relaxation rel2;
            THEN ("It can be move assigned with no leaks or exceptions") {
                REQUIRE_NOTHROW(rel = std::move(rel2));
            }
        }
    }
}

SCENARIO ("Black box testing of failures in constructing LP Relaxations",
          "[.LP][.Relaxaton][!shouldfail][valgrind]") {
    GIVEN ("The raw data structures in an LP Relaxation") {
        WHEN ("A nonzero rval occurs in the constructor") {
            THEN ("No memory is leaked") {
                int rval = 0;
        
                CPXLPptr cplex_lp = (CPXLPptr) NULL;
                CPXENVptr cplex_env = CPXopenCPLEX(&rval);

                rval = 5;

                if (rval) {
                    CPXcloseCPLEX(&cplex_env);
                    REQUIRE_NOTHROW(throw std::runtime_error("rval"));
                }
            }
        }
    }
}

#endif
