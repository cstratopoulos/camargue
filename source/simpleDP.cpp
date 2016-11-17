#include "simpleDP.hpp"
#include "DPgraph.hpp"
#include "PSEP_util.hpp"

#include <memory>

using std::vector;
using std::cout;
using std::cerr;

namespace PSEP {

Cut<dominoparity>::Cut(vector<int> &_delta,
		       vector<int> &_edge_marks,
		       vector<int> &_best_tour_nodes,
		       vector<int> &_perm,
		       SupportGraph &_G_s,
		       vector<int> &_support_elist,
		       vector<double> &_support_ecap) :
  candidates(_delta, _edge_marks,
	     _best_tour_nodes, _perm,
	     _G_s, _support_elist, _support_ecap)
{}

int Cut<dominoparity>::separate()
{
  int rval = 0;
  std::unique_ptr<DPCutGraph> witness;
  
  rval = candidates.get_light_teeth();
  if(rval) goto CLEANUP;

  candidates.weak_elim();
  

 CLEANUP:
  if(rval == 1)
    cerr << "Cut<dominoparity>::separate failed.\n";
  return rval;
}

}
