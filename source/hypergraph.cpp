#include "hypergraph.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>

using std::vector;

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



}
}
