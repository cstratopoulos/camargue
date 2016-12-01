#ifndef PSEP_CC_LPCUTS_HPP
#define PSEP_CC_LPCUTS_HPP

extern "C" {
#include <concorde/INCLUDE/tsp.h>
}

namespace PSEP {
namespace Cut {

class CCwrapper {
public:
  CCwrapper();
  ~CCwrapper();

  int cut_count() const { return _cutcount;}
  
  CCtsp_lpcut_in** pass_ptr() { return &cc_cut; }
  int* count_ptr() { return &_cutcount; }
  
private:
  CCtsp_lpcut_in *cc_cut;
  int _cutcount;
};

}
}

#endif
