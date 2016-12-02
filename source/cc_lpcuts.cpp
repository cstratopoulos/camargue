#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;

using lpcut_in = CCtsp_lpcut_in;

namespace PSEP {
namespace Cut {

CCwrapper::CCwrapper() : cc_cut(nullptr), cutcount(0) {}

CCwrapper::~CCwrapper()
{
  lpcut_in *it = cc_cut;
  
  while(it != nullptr){
    cc_cut = cc_cut->next;
    CCtsp_free_lpcut_in(it);
    delete(it);
    it = cc_cut;
  }
}

}
}
