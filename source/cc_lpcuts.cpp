#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;

using lpcut_in = CCtsp_lpcut_in;

namespace CMR {
namespace Cut {

LPcutIn::LPcutIn() : head_cut(nullptr), cutcount(0) {}

LPcutIn::~LPcutIn()
{
  lpcut_in *it = head_cut;
  
  while(it != nullptr){
    head_cut = head_cut->next;
    CCtsp_free_lpcut_in(it);
    delete(it);
    it = head_cut;
  }
}

void LPcutIn::filter_primal(CMR::TourGraph &TG)
{
  if(cutcount == 0 || begin() == nullptr) return;

  lpcut_in *current = begin();
  lpcut_in *prev = begin();

  while(current != end()){
    double cur_slack = CCtsp_cutprice(TG.pass_ptr(), current, TG.tour_array());

    if(cur_slack != 0){
      --cutcount;

      if(current == begin()){
	current = current->next;
	prev = current;
	CCtsp_free_lpcut_in(head_cut);
	delete(head_cut);
	head_cut = current;
      } else {
	prev->next = current->next;
	CCtsp_free_lpcut_in(current);
	delete(current);
	current = prev->next;
      }
      continue;
    }

    prev = current;
    current = current->next;
  }  
}

bool SegmentCuts::find_cuts()
{
  if(CCtsp_segment_cuts(cutq.pass_ptr(), cutq.count_ptr(),
			supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
			&supp_dat.support_elist[0],
			&supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_segment_cuts failed.");
  
  return(!cutq.empty());
}

bool BlockCombs::find_cuts()
{
  if(CCtsp_block_combs(cutq.pass_ptr(), cutq.count_ptr(),
		       supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
		       &supp_dat.support_elist[0],
		       &supp_dat.support_ecap[0], 1))
    throw std::runtime_error("CCtsp_block_combs failed.");

  if(cutq.empty())
    return false;

  cutq.filter_primal(TG);

  return(!cutq.empty());
}

bool FastBlossoms::find_cuts()
{
  if(CCtsp_fastblossom(cutq.pass_ptr(), cutq.count_ptr(),
		       supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
		       &supp_dat.support_elist[0],
		       &supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_fastblossom failed.");

  cutq.filter_primal(TG);
  if(!cutq.empty()) return true;

  if(CCtsp_ghfastblossom(cutq.pass_ptr(), cutq.count_ptr(),
			 supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
			 &supp_dat.support_elist[0],
			 &supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_ghfastblossom failed.");

  cutq.filter_primal(TG);
  return(!cutq.empty());
}

}
}
