#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;

using lpcut_in = CCtsp_lpcut_in;

namespace PSEP {
namespace Cut {

LPcutIn::LPcutIn() : cc_cut(nullptr), cutcount(0) {}

LPcutIn::~LPcutIn()
{
  lpcut_in *it = cc_cut;
  
  while(it != nullptr){
    cc_cut = cc_cut->next;
    CCtsp_free_lpcut_in(it);
    delete(it);
    it = cc_cut;
  }
}

void LPcutIn::filter_primal(PSEP::TourGraph &TG)
{
  if(cutcount == 0 || begin() == nullptr) return;
  
  lpcut_in *it = begin();

  while(it != nullptr){
    double slack = CCtsp_cutprice(TG.pass_ptr(), it, TG.tour_array());
    lpcut_in *del = it;
    it = it->next;

    if(slack != 0){
      --cutcount;
      
      if(begin() == del)
	cc_cut = begin()->next;

      if(del->next != nullptr)
	del->next->prev = del->prev;

      if(del->prev != nullptr)
	del->prev->next = del->next;

      CCtsp_free_lpcut_in(del);
      delete(del);
    }
  }
}

}
}
