#include "hypergraph.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>

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
    return pre_result;
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


}
}
