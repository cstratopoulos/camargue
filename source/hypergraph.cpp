#include "hypergraph.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

#include <cmath>

using std::unordered_map;
using std::vector;
using std::pair;

using std::cout;
using std::cerr;

using std::runtime_error;
using std::logic_error;
using std::exception;

using lpcut_in = CCtsp_lpcut_in;
using lpclique = CCtsp_lpclique;

namespace CMR {
namespace Sep {

HyperGraph::HyperGraph(CliqueBank &bank, const lpcut_in &cc_lpcut,
                       const vector<int> &tour) try :
    sense(cc_lpcut.sense), rhs(cc_lpcut.rhs), source_bank(&bank),
    source_toothbank(nullptr)
{
    for (int i = 0; i < cc_lpcut.cliquecount; ++i) {
        lpclique &cc_clq = cc_lpcut.cliques[i];
        cliques.push_back(source_bank->add_clique(cc_clq, tour));
    }
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("HyperGraph CC lpcut_in constructor failed.");
}

HyperGraph::HyperGraph(CliqueBank &bank, ToothBank &tbank,
                       const dominoparity &dp_cut, const double _rhs,
                       const std::vector<int> &tour) try :
    sense('L'), rhs(_rhs), source_bank(&bank), source_toothbank(&tbank)
{
    vector<int> nodes(dp_cut.degree_nodes);
    for(int &n : nodes)
        n = tour[n];
        
    cliques.push_back(source_bank->add_clique(nodes));

    for (const SimpleTooth &T : dp_cut.used_teeth)
        teeth.push_back(source_toothbank->add_tooth(T, tour));
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("HyperGraph dominoparity constructor failed.");
}

HyperGraph::~HyperGraph()
{
    if (source_bank != nullptr)
        for (Clique::Ptr &ref : cliques)
            source_bank->del_clique(ref);

    if (source_toothbank != nullptr)
        for (Tooth::Ptr &ref : teeth)
            source_toothbank->del_tooth(ref);
}

double HyperGraph::get_coeff(int end0, int end1) const
{
    if (end0 == end1)
        throw logic_error("Edge has same endpoints in HyperGraph::get_coeff.");
    
    double result = 0.0;
    
    if (cut_type() == Type::Standard) {
        const vector<int> &perm = source_bank->ref_perm();

        int end0_ind = perm[end0];
        int end1_ind = perm[end1];

        for (const Clique::Ptr &clq_ref : cliques) {
            bool contains_end0 = clq_ref->contains(end0_ind);
            bool contains_end1 = clq_ref->contains(end1_ind);

            if (contains_end0 != contains_end1)
                result += 1.0;
        }

        return result;
    }

    //else it is a Simple DP
    int pre_result = 0;
    //get handle coeffs
    const vector<int> &handle_perm = source_bank->ref_perm();
    int end0_ind = handle_perm[end0];
    int end1_ind = handle_perm[end1];
    
    const Clique::Ptr &handle_clq = cliques[0];
    
    bool contains_end0 = handle_clq->contains(end0_ind);
    bool contains_end1 = handle_clq->contains(end1_ind);
    
    if (contains_end0 && contains_end1) //in E(H)
        pre_result += 2;
    else if (contains_end0 != contains_end1) // in delta(H)
        pre_result += 1;

    const vector<int> &tooth_perm = source_toothbank->ref_perm();
    
    end0_ind = tooth_perm[end0];
    end1_ind = tooth_perm[end1];
    contains_end0 = false;
    contains_end1 = false;

    for (const Tooth::Ptr &tooth : teeth) {
        const Clique &root_clq = tooth->set_pair()[0];
        const Clique &bod_clq = tooth->set_pair()[1];

        bool root_end0 = false;
        bool root_end1 = false;

        if ((root_end0 = root_clq.contains(end0_ind)))
            if (bod_clq.contains(end1_ind)) {
                pre_result += 1;
                continue;
            }

        if ((root_end1 = root_clq.contains(end1_ind)))
            if (bod_clq.contains(end0_ind)) {
                pre_result += 1;
                continue;
            }

        if (root_end0 || root_end1)
            continue;

        if (bod_clq.contains(end0_ind) && bod_clq.contains(end1_ind))
            pre_result += 2;
    }

    pre_result /= 2;
    return static_cast<double>(pre_result);
}

void HyperGraph::get_coeffs(const std::vector<Price::PrEdge> &edges,
                            vector<int> &rmatind,
                            vector<double> &rmatval) const
{
    rmatind.clear();
    rmatval.clear();

    const vector<int> &def_tour = source_bank->ref_tour();
    int ncount = def_tour.size();
    
    std::map<int, double> coeff_map;
    vector<bool> node_marks(ncount, false);

    if (cut_type() == Type::Standard) {
        for (const Clique::Ptr &clq_ref : cliques) {
            for (const Segment &seg : clq_ref->seg_list())
                for (int k = seg.start; k <= seg.end; ++k)
                    node_marks[def_tour[k]] = true;

            for (int i = 0; i < edges.size(); ++i) {
                if (node_marks[edges[i].end[0]] !=
                    node_marks[edges[i].end[1]]) {
                    if (coeff_map.count(i))
                        coeff_map[i] += 1.0;
                    else
                        coeff_map[i] = 1.0;
                }
            }

            node_marks = vector<bool>(ncount, false);
        }

        rmatind.reserve(coeff_map.size());
        rmatval.reserve(coeff_map.size());

        for (std::pair<const int, double> &kv : coeff_map) {
            rmatind.push_back(kv.first);
            rmatval.push_back(kv.second);
        }

        return;
    } //else it is a domino cut

    const Clique::Ptr &handle_ref = cliques[0];
    for (const Segment &seg : handle_ref->seg_list())
        for (int k = seg.start; k <= seg.end; ++k)
            node_marks[def_tour[k]] = true;

    for (int i = 0; i < edges.size(); ++i) {
        bool e0 = node_marks[edges[i].end[0]];
        bool e1 = node_marks[edges[i].end[1]];

        if (e0 && e1) {
            if (coeff_map.count(i)) coeff_map[i] += 2.0;
            else coeff_map[i] = 2.0;
        } else if (e0 != e1) {
            if (coeff_map.count(i)) coeff_map[i] += 1.0;
            else coeff_map[i] = 1.0;
        }
    }

    node_marks = vector<bool>(ncount, false);

    const vector<int> &tooth_perm = source_toothbank->ref_perm();

    for (const Tooth::Ptr &t_ref : teeth) {
        const Clique &root_clq = t_ref->set_pair()[0];
        const Clique &bod_clq = t_ref->set_pair()[1];

        for (const Segment &seg : bod_clq.seg_list()) 
            for (int k = seg.start; k <= seg.end; ++k)
                node_marks[def_tour[k]] = true;

        for (int i = 0; i < edges.size(); ++i) {
            int e0 = edges[i].end[0];
            int e1 = edges[i].end[1];
            
            bool bod_e0 = node_marks[e0];
            bool bod_e1 = node_marks[e1];

            if (bod_e0 && bod_e1) {
                if (coeff_map.count(i)) coeff_map[i] += 2.0;
                else coeff_map[i] = 2.0;
            } else if (bod_e0 != bod_e1) {
                if (bod_e0) {
                    if (root_clq.contains(tooth_perm[e1])) {
                        if (coeff_map.count(i)) coeff_map[i] += 1.0;
                        else coeff_map[i] = 1.0;
                    }
                } else { //bod_e1
                    if (root_clq.contains(tooth_perm[e0])) {
                        if (coeff_map.count(i)) coeff_map[i] += 1.0;
                        else coeff_map[i] = 1.0;
                    }
                }
            }
        }

        node_marks = vector<bool>(ncount, false);
    }

    rmatind.reserve(coeff_map.size());
    rmatval.reserve(coeff_map.size());

    for (std::pair<const int, double> &kv : coeff_map) {
        rmatind.push_back(kv.first);
        rmatval.push_back(kv.second);
    }

    for (double &coeff : rmatval)
        if (fabs(coeff >= Epsilon::Zero)) {
            coeff /= 2;
            coeff = floor(coeff);
        }

    
}

ExternalCuts::ExternalCuts(const vector<int> &tour, const vector<int> &perm)
try : node_count(tour.size()), clique_bank(tour, perm), tooth_bank(tour, perm)
{
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("ExternalCuts constructor failed.");
}

void ExternalCuts::add_cut(const lpcut_in &cc_lpcut,
                           const vector<int> &current_tour)
{
    cuts.emplace_back(clique_bank, cc_lpcut, current_tour);
}

void ExternalCuts::add_cut(const dominoparity &dp_cut, const double rhs,
                           const vector<int> &current_tour)
{
    cuts.emplace_back(clique_bank, tooth_bank, dp_cut, rhs, current_tour);
}

void ExternalCuts::del_cuts(const vector<int> &delset)
{
    int i = 0;

    for (HyperGraph &H : cuts) {
        if (delset[i + node_count] == -1)
            H.rhs = '\0';
        ++i;
    }

    cuts.erase(std::remove_if(cuts.begin(), cuts.end(),
                              [](const HyperGraph &H) {
                                  return H.rhs == '\0';
                              }),
               cuts.end());
}

void ExternalCuts::get_col(const int end0, const int end1,
                           vector<int> &cmatind, vector<double> &cmatval)
{
    runtime_error err("Problem in ExternalCuts::get_col");
    
    if (end0 == end1) {
        cerr << "Edge has same endpoints.\n";
        throw err;
    }
    
    cmatind.clear();
    cmatval.clear();
    
    int lp_size = node_count + cuts.size();

    try {
        cmatind.reserve(lp_size);
        cmatval.reserve(lp_size);

        cmatind.push_back(end0);
        cmatind.push_back(end1);

        cmatval.push_back(1.0);
        cmatval.push_back(1.0);

        for (int i = 0; i < cuts.size(); ++i) {
            int index = i + node_count;
            double coeff = cuts[i].get_coeff(end0, end1);

            if (coeff != 0.0) {
                cmatind.push_back(index);
                cmatval.push_back(coeff);
            }
        }
    } CMR_CATCH_PRINT_THROW("Couldn't push back column coeffs/inds", err);
}

void ExternalCuts::get_duals(const LP::Relaxation &relax,
                             vector<double> &node_pi,
                             vector<double> &node_pi_est,
                             vector<double> &cut_pi,
                             unordered_map<Clique, double> clique_pi) const
{
    runtime_error err("Problem in ExternalCuts::get_duals.");
    using CutType = HyperGraph::Type;

    if (relax.num_rows() != node_count + cuts.size()) {
        cerr << "Relaxation row count: " << relax.num_rows() << ", "
             << "ExternalCuts expects: " << (node_count + cuts.size()) << "\n";
        throw logic_error("Size mismatch in ExternalCuts::get_duals.");
    }

    clique_pi.clear();

    try {
        clique_pi.reserve(clique_bank.size());

        relax.get_pi(node_pi, 0, node_count - 1);
        node_pi_est = node_pi;

        if (!cuts.empty())
            relax.get_pi(cut_pi, node_count, relax.num_rows() - 1);
    } CMR_CATCH_PRINT_THROW("Couldn't get/allocate pi containers.", err);

    //get clique_pi for non-domino cuts
    for (int i = 0; i < cuts.size(); ++i) {
        const Sep::HyperGraph &H = cuts[i];
        if (H.cut_type() != CutType::Standard)
            continue;

        double pival = cut_pi[i];

        for (const Sep::Clique::Ptr &clq_ref : H.get_cliques()) {
            if (clique_pi.count(*clq_ref) == 0)
                clique_pi[*clq_ref] = 0.0;
            
            clique_pi[*clq_ref] += pival;
        }
    }

    const std::vector<int> &def_tour = clique_bank.ref_tour();

    //use clique_pi to build node_pi for all cliques with nonzero pi
    for (const std::pair<Sep::Clique, double> &kv : clique_pi) {
        const Sep::Clique &clq = kv.first;
        double pival = kv.second;

        if (pival > 0.0) {
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];
                    
                    node_pi[node] += pival;
                    node_pi_est[node] += pival;
                }
            }            
        } else if (pival < 0.0) {
            for (const Segment &seg : clq.seg_list()) {
                for (int k = seg.start; k <= seg.end; ++k) {
                    int node = def_tour[k];
                    
                    node_pi[node] += pival;
                }
            }
        }
    }

    //now get node_pi_est for domino cuts, skipping standard ones
    for (int i = 0; i < cuts.size(); ++i) {
        const Sep::HyperGraph &H = cuts[i];
        if (H.cut_type() != CutType::Domino)
            continue;

        double pival = cut_pi[i];
        if (pival <= 0.0)
            continue;
        
        const Sep::Clique::Ptr &handle_ref = H.get_cliques()[0];
        for (const Segment &seg : handle_ref->seg_list()) {
            for (int k = seg.start; k <= seg.end; ++k) {
                int node = def_tour[k];
                node_pi_est[node] += pival;
            }
        }

        for (const Sep::Tooth::Ptr &T : H.get_teeth())
            for (const Sep::Clique &tpart : T->set_pair())
                for (const Segment &seg : tpart.seg_list())
                    for (int k = seg.start; k <= seg.end; ++k) {
                        int node = def_tour[k];

                        node_pi_est[node] += pival;
                    }
    }    
    
}


}
}
