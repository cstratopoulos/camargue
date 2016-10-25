#include "DPgraph.hpp"

#include <iostream>

using std::vector;

namespace PSEP {

DPCutGraph::DPCutGraph(const vector<vector<SimpleTooth::Ptr>> &_teeth) :
  light_teeth(_teeth) { CCcut_GHtreeinit(&gh_tree); }

DPCutGraph::~DPCutGraph(){ CCcut_GHtreefree(&gh_tree); }

}
