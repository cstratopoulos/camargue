#include "tests.hpp"

#ifdef CMR_DO_TESTS

#include "hypergraph.hpp"
#include "datagroups.hpp"
#include "separator.hpp"
#include "util.hpp"

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <iostream>

#include <catch.hpp>

using std::min;
using std::max;
using std::abs;
using std::array;
using std::vector;

using std::string;
using std::cout;

SCENARIO ("Comparing HyperGraph edge coeffs to sparse rows",
          "[Sep][HyperGraph][get_coeff][get_coeffs]") {
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

        CMR::Data::GraphGroup g_dat;
        CMR::Data::BestGroup b_dat;
        CMR::Data::SupportGroup s_dat;
        vector<double> lp_edges;
        CMR::Data::Instance inst;
        CMR::Data::KarpPartition kpart;

        GIVEN ("A subtour lp solution for " + fname) {
            CMR::Data::make_cut_test(probfile, solfile, subtourfile, g_dat,
                                     b_dat, lp_edges, s_dat, inst);
            int ncount = g_dat.core_graph.node_count();

            
            CMR::Sep::CutTranslate translator(g_dat);
            
            kpart = CMR::Data::KarpPartition(inst);

            CMR::Sep::Separator sep(g_dat, b_dat, s_dat, kpart,
                                    CMR::IntMax);
            CMR::TourGraph TG(b_dat.best_tour_edges,
                              g_dat.core_graph.get_edges(),
                              b_dat.perm);

            CMR::Sep::CliqueBank cbank(b_dat.best_tour_nodes,
                                       b_dat.perm);
            CMR::Sep::ToothBank tbank(cbank);

            vector<CMR::Price::edge> pr_edges(g_dat.core_graph.edge_count());
            for (int i = 0; i < pr_edges.size(); ++i)
                pr_edges[i].end = g_dat.core_graph.get_edge(i).end;

            WHEN ("Cuts are found") {
                REQUIRE(sep.find_cuts(TG));

                THEN ("All the sparse row coeffs match HyperGraph coeffs.") {
                    std::array<CMR::Sep::LPcutList*,
                               2> qlist{&sep.fast2m_q,
                                        &sep.blkcomb_q};
                    
                    vector<int> &perm = b_dat.perm;
                    vector<int> &tour = b_dat.best_tour_nodes;

                    for (CMR::Sep::LPcutList *qptr : qlist) {
                        if (qptr->empty())
                            continue;

                        for (CCtsp_lpcut_in *cur = qptr->begin(); cur;
                             cur = cur->next) {
                            vector<int> rmatind;
                            vector<double> rmatval;
                            char sense;
                            double rhs;

                            translator.get_sparse_row(*cur, perm, rmatind,
                                                      rmatval, sense, rhs);

                            CMR::Sep::HyperGraph hg(cbank, *cur, tour);

                            for (int i = 0; i < rmatind.size(); ++i) {
                                int edge_ind = rmatind[i];
                                double ref_coeff = rmatval[i];
                                CMR::Edge e =
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

                    if (!sep.dp_q.empty()) {
                        cout << "\tFOUND DP CUTS\n";
                        for (auto it = sep.dp_q.begin();
                             it != sep.dp_q.end(); ++it) {
                            CMR::Sep::dominoparity &dp_cut = *it;
                            vector<int> rmatind;
                            vector<double> rmatval;
                            char sense;
                            double rhs;

                            translator.get_sparse_row(dp_cut, tour, rmatind,
                                                      rmatval, sense, rhs);

                            CMR::Sep::HyperGraph hg(cbank, tbank, dp_cut, rhs,
                                                    tour);
                            
                            for (int i = 0; i < rmatind.size(); ++i) {
                                int edge_ind = rmatind[i];
                                double ref_coeff = rmatval[i];
                                CMR::Edge e =
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
