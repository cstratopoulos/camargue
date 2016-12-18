#include "hypergraph.hpp"

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
    sense(cc_lpcut.sense), rhs(cc_lpcut.rhs), source_bank(bank)
{
    for (int i = 0; i < cc_lpcut.cliquecount; ++i) {
        lpclique &cc_clq = cc_lpcut.cliques[i];
        cliques.push_back(source_bank.add_clique(cc_clq, tour));
    }
    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Hypergraph CC lpcut_in constructor failed.");
}

HyperGraph::~HyperGraph()
{
    for (Clique::Ptr &ref : cliques)
        source_bank.del_clique(ref);
}

/*
DominoCut::DominoCut(CliqueBank &bank, CMR::Sep::dominoparity &dp_cut,
                     const vector<int> &tour) try :
    source_bank(bank)
{
    for (int &i : dp_cut.degree_nodes)
        i = tour[i];

    for (auto &ends : dp_cut.nonneg_edges) {
        ends.first = tour[ends.first];
        ends.second = tour[ends.second];
    }

    nonneg_edges = dp_cut.nonneg_edges;
    
    for (CMR::SimpleTooth &T : dp_cut.used_teeth) {
        teeth.push_back(pair<int, Clique::Ptr>(tour[T.root],
                                               bank.add_clique(T.body_start,
                                                               T.body_end,
                                                               tour)));
    }
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Problem in DominoCut constructor.");
}


DominoCut::~DominoCut()
{
    source_bank.del_clique(handle);
    for(auto &tooth : teeth)
        source_bank.del_clique(tooth.second);
}
*/

ExternalCuts::ExternalCuts(const vector<int> &tour, const vector<int> &perm)
try : next_row(tour.size()), clique_bank(tour, perm) {
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("ExternalCuts constructor failed.");
}

void ExternalCuts::add_cut(const lpcut_in &cc_lpcut,
                           const vector<int> &current_tour)
{
    cuts.insert(std::pair<int, HyperGraph>(next_row++,
					   HyperGraph(clique_bank, cc_lpcut,
						      current_tour)));
}



}
}
