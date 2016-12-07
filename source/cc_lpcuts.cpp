#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>
#include <utility>

using std::cout;
using std::endl;
using std::cerr;

using lpcut_in = CCtsp_lpcut_in;

namespace CMR {
namespace Cut {

LPcutList::LPcutList() noexcept : head_cut(), cutcount(0) {}

LPcutList::LPcutList(CCtsp_lpcut_in *head, int count) noexcept :
  head_cut(head), cutcount(count) {}

LPcutList::LPcutList(LPcutList &&L) noexcept :
  head_cut(std::move(L.head_cut)), cutcount(L.cutcount) { L.cutcount = 0; }

LPcutList &LPcutList::operator=(LPcutList &&L) noexcept {
  head_cut = std::move(L.head_cut);
  cutcount = L.cutcount;
  L.cutcount = 0;
  return *this;
}

void LPcutList::filter_primal(CMR::TourGraph &TG)
{
  if(cutcount == 0 || !head_cut) return;

  lpcut_in *current = head_cut.get();
  lpcut_in *prev = current;

  while(current){
    double slack = CCtsp_cutprice(TG.pass_ptr(), current, TG.tour_array());

    if(slack != 0){
      --cutcount;

      if(current == head_cut.get()){
	current = head_cut->next;
	head_cut->next = nullptr;
	head_cut.reset(current);
	current = head_cut.get();
	prev = current;
      } else {
	prev->next = current->next;
	CCtsp_free_lpcut_in(current);
	CC_IFFREE(current, lpcut_in);
	current = prev->next;
      }
    } else {
      prev = current;
      current = current->next;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//                        CONCORDE SEPARATORS                                //
///////////////////////////////////////////////////////////////////////////////

bool SegmentCuts::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;
  
  if(CCtsp_segment_cuts(&head, &cutcount,
			supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
			&supp_dat.support_elist[0],
			&supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_segment_cuts failed.");

  if(cutcount == 0) return false;

  cutq = LPcutList(head, cutcount);
  cutq.filter_primal(TG);
  
  return(!cutq.empty());
}

bool BlockCombs::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;
  
  if(CCtsp_block_combs(&head, &cutcount,
		       supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
		       &supp_dat.support_elist[0],
		       &supp_dat.support_ecap[0], 1))
    throw std::runtime_error("CCtsp_block_combs failed.");

  if(cutcount == 0) return false;

  cutq = LPcutList(head, cutcount);
  cutq.filter_primal(TG);

  return(!cutq.empty());
}

bool FastBlossoms::find_cuts()
{
  int cutcount = 0;
  lpcut_in *head = (lpcut_in *) NULL;
  
  if(CCtsp_fastblossom(&head, &cutcount,
		       supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
		       &supp_dat.support_elist[0],
		       &supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_fastblossom failed.");

  if(cutcount == 0) return false;

  cutq = LPcutList(head, cutcount);
  cutq.filter_primal(TG);
  if(!cutq.empty()) return true;

  cutcount = 0;
  head = (lpcut_in *) NULL;

  if(CCtsp_ghfastblossom(&head, &cutcount,
			 supp_dat.G_s.node_count, supp_dat.G_s.edge_count,
			 &supp_dat.support_elist[0],
			 &supp_dat.support_ecap[0]))
    throw std::runtime_error("CCtsp_ghfastblossom failed.");

  cutq = LPcutList(head, cutcount);
  cutq.filter_primal(TG);
  
  return(!cutq.empty());
}

}
}
