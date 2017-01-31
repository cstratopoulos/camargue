#include "config.hpp"

#ifdef CMR_DO_TESTS

#include "hypergraph.hpp"
#include "datagroups.hpp"
#include "separator.hpp"
#include "solver.hpp"
#include "util.hpp"

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <functional>

#include <catch.hpp>

using std::min;
using std::max;
using std::abs;
using std::array;
using std::vector;

using std::string;
using std::cout;

SCENARIO ("Comparing HyperGraph edge coeffs to CPLEX coefs",
          "[HyperGraph][get_coeff][Sep][LP]")
{
    using namespace CMR;
    vector<string> probs {
        "ulysses16",
        "dantzig42",
        "eil51",
        "rat99",
        "lin318",
        "d493",
        "p654",
        "u724",
        "dsj1000",
        "pr1002",
        "d2103",
        };

    for (string &fname : probs) {
        GIVEN ("A priceless cutting_loop run on " + fname) {
            WHEN ("We get coefficients of individual edges") {
                THEN ("They agree with those from CPLEX") {
                    string probfile = "problems/" + fname + ".tsp";

                    OutPrefs outprefs;        
                    Solver solver(probfile, 99, outprefs);
                    
                    solver.cutting_loop(false);

                    const Data::GraphGroup &g_dat = solver.graph_info();
                    const LP::CoreLP &core_lp = solver.get_core_lp();

                    const vector<Graph::Edge> &edges = g_dat.core_graph
                    .get_edges();
                    
                    int ncount = g_dat.core_graph.node_count();
                    int numrows = core_lp.num_rows();

                    for (int i = 0; i < edges.size(); ++i) {
                        vector<int> ex_cmatind;
                        vector<double> ex_cmatval;
                        vector<int> cpx_cmatind;
                        vector<double> cpx_cmatval;

                        REQUIRE_NOTHROW(core_lp.external_cuts()
                                        .get_col(edges[i].end[0],
                                                 edges[i].end[1],
                                                 ex_cmatind,
                                                 ex_cmatval));
                        REQUIRE_NOTHROW(core_lp.get_col(i, cpx_cmatind,
                                                        cpx_cmatval));

                        REQUIRE(ex_cmatind.size() == cpx_cmatind.size());
                        REQUIRE(ex_cmatind == cpx_cmatind);
                        REQUIRE(ex_cmatval == cpx_cmatval);
                    }
                }
            }
        }
    }
}


SCENARIO ("Comparing HyperGraph edge coeffs to comb/domino sparse rows",
          "[Sep][HyperGraph][get_coeff][get_coeffs]") {
    using namespace CMR;
    vector<string> probs{
        "dantzig42",
        "pr76",
        "lin105",
        "d493",
    };

    for (string &fname : probs) {
        string
        probfile = "problems/" + fname + ".tsp",
        solfile = "test_data/tours/" + fname + ".sol",
        subtourfile = "test_data/subtour_lp/" + fname + ".sub.x";

        Data::GraphGroup g_dat;
        Data::BestGroup b_dat;
        Data::SupportGroup s_dat;
        vector<double> lp_edges;
        Data::Instance inst;
        Data::KarpPartition kpart;

        GIVEN ("A subtour lp solution for " + fname) {
            Data::make_cut_test(probfile, solfile, subtourfile, g_dat,
                                     b_dat, lp_edges, s_dat, inst);
            int ncount = g_dat.core_graph.node_count();

            
            Sep::CutTranslate translator(g_dat);
            
            kpart = Data::KarpPartition(inst);
            Graph::TourGraph TG(b_dat.best_tour_edges,
                                     g_dat.core_graph.get_edges(),
                                     b_dat.perm);

            Sep::Separator sep(g_dat, b_dat, s_dat, kpart, TG,
                                    IntMax);


            Sep::CliqueBank cbank(b_dat.best_tour_nodes,
                                       b_dat.perm);
            Sep::ToothBank tbank(cbank);

            vector<Price::PrEdge<double>> pr_edges(g_dat.core_graph.edge_count());
            for (int i = 0; i < pr_edges.size(); ++i)
                pr_edges[i].end = g_dat.core_graph.get_edge(i).end;

            WHEN ("Cuts are found") {
                bool fast2m = sep.fast2m_sep();
                bool blkcomb = sep.blkcomb_sep();
                bool dp = sep.simpleDP_sep();
                bool found = fast2m || blkcomb || dp;
                REQUIRE(found);

                THEN ("All the sparse row coeffs match HyperGraph coeffs.") {
                    std::array<const CCtsp_lpcut_in*, 2>
                    qheads{sep.fastblossom_q().begin(),
                        sep.blockcomb_q().begin()};
                                        
                    vector<int> &perm = b_dat.perm;
                    vector<int> &tour = b_dat.best_tour_nodes;

                    for (const CCtsp_lpcut_in* headptr : qheads)
                        for (const CCtsp_lpcut_in *cur = headptr;
                             cur; cur = cur->next) {
                            vector<int> rmatind;
                            vector<double> rmatval;
                            char sense;
                            double rhs;

                            translator.get_sparse_row(*cur, perm, rmatind,
                                                      rmatval, sense, rhs);

                            Sep::HyperGraph hg(cbank, *cur, tour);

                            for (int i = 0; i < rmatind.size(); ++i) {
                                int edge_ind = rmatind[i];
                                double ref_coeff = rmatval[i];
                                Graph::Edge e =
                                g_dat.core_graph.get_edge(edge_ind);

                                double hg_coeff = hg.get_coeff(e.end[0],
                                                               e.end[1]);

                                REQUIRE(ref_coeff == hg_coeff);
                            }

                            vector<int> pr_rmatind;
                            vector<double> pr_rmatval;
                            hg.get_coeffs(pr_edges , pr_rmatind, pr_rmatval);
                            REQUIRE(pr_rmatind == rmatind);
                            REQUIRE(pr_rmatval == rmatval);
                        }
                    

                    if (!sep.simpleDP_q().empty()) {
                        cout << "\tFOUND DP CUTS\n";
                        for (auto it = sep.simpleDP_q().begin();
                             it != sep.simpleDP_q().end(); ++it) {
                            const Sep::dominoparity &dp_cut = *it;
                            vector<int> rmatind;
                            vector<double> rmatval;
                            char sense;
                            double rhs;

                            translator.get_sparse_row(dp_cut, tour, rmatind,
                                                      rmatval, sense, rhs);

                            Sep::HyperGraph hg(cbank, tbank, dp_cut, rhs,
                                                    tour);
                            
                            for (int i = 0; i < rmatind.size(); ++i) {
                                int edge_ind = rmatind[i];
                                double ref_coeff = rmatval[i];
                                Graph::Edge e =
                                g_dat.core_graph.get_edge(edge_ind);

                                double hg_coeff = hg.get_coeff(e.end[0],
                                                               e.end[1]);

                                REQUIRE(ref_coeff == hg_coeff);
                            }

                            vector<int> pr_rmatind;
                            vector<double> pr_rmatval;
                            hg.get_coeffs(pr_edges , pr_rmatind, pr_rmatval);
                            REQUIRE(pr_rmatind == rmatind);
                            REQUIRE(pr_rmatval == rmatval);
                        }
                    }
                }
            }

        }
    }
}

#endif //CMR_DO_TESTS
