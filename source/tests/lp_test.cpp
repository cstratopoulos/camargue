#include "config.hpp"

#ifdef CMR_DO_TESTS

#include <cplex.h>

#include "fixed64.hpp"
#include "lp_interface.hpp"
#include "core_lp.hpp"
#include "separator.hpp"

#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <catch.hpp>

using std::abs;
using std::cout;
using std::vector;
using std::unique_ptr;
using std::string;
using std::pair;

SCENARIO ("Benchmarking rounds of cuts",
          "[LP][CoreLP][primal_pivot][Sep][Separator][benchmark]") {
    using namespace CMR;
    vector<string> probs{"pr76", "a280", "p654", "pr1002", "rl1304"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            THEN ("We test the changes from adding cuts in rounds") {
                    Data::Instance inst("problems/" + prob + ".tsp", 99);
                    Data::GraphGroup g_dat(inst);
                    Data::BestGroup b_dat(inst, g_dat);
                    LP::CoreLP core(g_dat, b_dat);

                    vector<double> tourx = core.lp_vec();
                    double tourlen = core.get_objval();

                    REQUIRE_NOTHROW(core.primal_pivot());

                    vector<double> pivx = core.lp_vec();
                    double pval = core.get_objval();
                    cout << "\tFirst pivot val: " << pval << "\n";

                    Data::SupportGroup s_dat;
                    int ncount = inst.node_count();

                    s_dat.reset(ncount,
                                g_dat.core_graph.get_edges(), pivx,
                                g_dat.island);

                    Graph::TourGraph TG(b_dat.best_tour_edges,
                                        g_dat.core_graph.
                                        get_edges(), b_dat.perm);
                    Data::KarpPartition kpart;

                    unique_ptr<Sep::Separator> sep =
                    util::make_unique<Sep::Separator>(g_dat, b_dat, s_dat,
                                                      kpart, TG);

                    if (sep->connect_sep()) {
                        core.pivot_back();
                        core.add_cuts(sep->connect_cuts_q());
                        LP::PivType piv = core.primal_pivot();
                        vector<double> newx = core.lp_vec();
                        double newpiv = core.get_objval();
                        double delta = abs(newpiv - pval);
                        cout << "\tConnect pivot val: " << newpiv << "\n"
                             << "\t\tDelta " << delta << "\t"
                             << (delta / tourlen) << " tour ratio\n"
                             << "\tStat " << piv << "\n";
                        pivx = newx;
                        pval = newpiv;
                        s_dat.reset(ncount, g_dat.core_graph.get_edges(), pivx,
                                    g_dat.island);
                        sep = (util::make_unique<Sep::Separator>(g_dat,
                                                                    b_dat,
                                                                    s_dat,
                                                                    kpart,
                                                                    TG));
                    }

                    if (sep->fast2m_sep()) {
                        core.pivot_back();
                        core.add_cuts(sep->fastblossom_q());
                        LP::PivType piv = core.primal_pivot();
                        vector<double> newx = core.lp_vec();
                        double newpiv = core.get_objval();
                        double delta = abs(newpiv - pval);
                        cout << "\tIs fast2m piv dif from first piv: "
                             << (newx == pivx) << "\n";
                        cout << "\tFast2m pivot val: " << newpiv << "\n"
                             << "\t\tDelta " << delta << "\t"
                             << (delta / tourlen) << " tour ratio\n"
                             << "\tStat " << piv << "\n";
                        pivx = newx;
                        pval = newpiv;
                        s_dat.reset(ncount, g_dat.core_graph.get_edges(), pivx,
                                    g_dat.island);
                        sep = (util::make_unique<Sep::Separator>(g_dat,
                                                                    b_dat,
                                                                    s_dat,
                                                                    kpart,
                                                                    TG));
                    }

                    if (sep->blkcomb_sep()) {
                        core.pivot_back();
                        core.add_cuts(sep->blockcomb_q());
                        LP::PivType piv = core.primal_pivot();
                        vector<double> newx = core.lp_vec();
                        double newpiv = core.get_objval();
                        double delta = abs(newpiv - pval);
                        cout << "\tBlkcomb pivot val: " << newpiv << "\n"
                             << "\t\tDelta " << delta << "\t"
                             << (delta / tourlen) << " tour ratio\n"
                             << "\tStat " << piv << "\n";
                        pivx = newx;
                        pval = newpiv;
                        s_dat.reset(ncount, g_dat.core_graph.get_edges(), pivx,
                                    g_dat.island);
                        sep = (util::make_unique<Sep::Separator>(g_dat,
                                                                    b_dat,
                                                                    s_dat,
                                                                    kpart,
                                                                    TG));
                    }
            }
        }
    }
}

SCENARIO ("Performing single pivots",
          "[LP][CoreLP][primal_pivot][Sep][Separator]") {
    using namespace CMR;
    vector<string> probs{"pr76", "a280", "p654", "pr1002", "rl1304"};

    for (string &prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("We pivot to a new solution") {
                THEN ("The LP vector changes") {
                    Data::Instance inst("problems/" + prob + ".tsp", 99);
                    Data::GraphGroup g_dat(inst);
                    Data::BestGroup b_dat(inst, g_dat);
                    LP::CoreLP core(g_dat, b_dat);

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

                            Data::SupportGroup s_dat;

                            s_dat.reset(inst.node_count(),
                                        g_dat.core_graph.get_edges(), pivx,
                                        g_dat.island);

                            Graph::TourGraph TG(b_dat.best_tour_edges,
                                                     g_dat.core_graph.
                                                     get_edges(), b_dat.perm);
                            Data::KarpPartition kpart;
                            Sep::Separator sep(g_dat, b_dat, s_dat,
                                                    kpart, TG);

                            bool found_some = false;

                            if (sep.segment_sep()) {
                                found_some = true;
                                core.add_cuts(sep.segment_q());
                            } else if (sep.fast2m_sep()) {
                                found_some = true;
                                core.add_cuts(sep.fastblossom_q());
                            } else if (sep.blkcomb_sep()) {
                                found_some = true;
                                core.add_cuts(sep.blockcomb_q());
                            } else if (sep.connect_sep()) {
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
    using namespace CMR;
    vector<string> probs{"dantzig42", "lin318", "d493", "pr1002", "pcb3038"};

    for (string& prob : probs) {
        GIVEN ("The TSP instance " + prob) {
            WHEN ("A core LP is constructed") {
                Data::Instance inst("problems/" + prob + ".tsp", 99);
                Data::GraphGroup g_dat(inst);
                Data::BestGroup b_dat(inst, g_dat);
                THEN ("Its constructor doesn't throw.") {
                    std::unique_ptr<LP::CoreLP> core;
                    REQUIRE_NOTHROW(core =
                                    util::make_unique<LP::CoreLP>(g_dat,
                                                                  b_dat));
                }

                AND_THEN ("The degree LP is feasible at the best tour") {
                    LP::CoreLP core(g_dat, b_dat);
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
    using namespace CMR;
    WHEN("A relaxation is constructed"){
        THEN("Its constructor doesn't throw"){
            std::unique_ptr<LP::Relaxation> rel;
            REQUIRE_NOTHROW(rel = util::make_unique<LP::Relaxation>());
        }
    }

    GIVEN ("A constructed Relaxation") {
        LP::Relaxation rel;
        WHEN ("Another one is constructed") {
            LP::Relaxation rel2;
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
