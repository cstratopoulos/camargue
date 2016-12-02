#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;

using lpcut_in = CCtsp_lpcut_in;

namespace PSEP {
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

void LPcutIn::filter_primal(PSEP::TourGraph &TG)
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

}
}
