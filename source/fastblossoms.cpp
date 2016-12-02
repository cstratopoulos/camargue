#include "fastblossoms.hpp"

#include <stdexcept>

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

namespace PSEP {
namespace Cut {

FastBlossoms::FastBlossoms(PSEP::Data::GraphGroup &_graph_dat,
			   PSEP::Data::BestGroup &_best_dat,
			   PSEP::Data::SupportGroup &_supp_dat,
			   PSEP::TourGraph &_TG,
			   LPcutIn &_cutq) :
  graph_dat(_graph_dat),
  best_dat(_best_dat),
  supp_dat(_supp_dat),
  TG(_TG), cutq(_cutq) {}

bool FastBlossoms::find_cuts()
{
  if(find_odd())
    return true;

  return find_gh();
}

bool FastBlossoms::find_odd()
{
  if(CCtsp_fastblossom(cutq.pass_ptr(), cutq.count_ptr(),
		       supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
		       &supp_dat.support_elist[0],
		       &supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_fastblossom failed.");

  if(cutq.empty()){ return false; }

  // for(auto it = cutq.begin(); it != cutq.end(); it = it->next){
  //   double slack = CCtsp_cutprice(TG.pass_ptr(), it, TG.tour_array());
    
  //   if(slack != 0) cout << "\tvvvv NON TIGHT CUT vvvvvv\n";
  //   CCtsp_print_lpcut_in(it);
  //   it = it->next;
  // }
  
  cutq.filter_primal(TG);

  return(!cutq.empty());
}

bool FastBlossoms::find_gh()
{
  if(CCtsp_ghfastblossom(cutq.pass_ptr(), cutq.count_ptr(),
		       supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
		       &supp_dat.support_elist[0],
		       &supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_ghfastblossom failed.");

  if(cutq.empty())
    return false;

  cutq.filter_primal(TG);
  
  return(!cutq.empty());
}

}
}
