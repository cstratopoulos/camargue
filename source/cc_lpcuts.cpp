#include "cc_lpcuts.hpp"

#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;

using lpcut_in = CCtsp_lpcut_in;

namespace PSEP {
namespace Cut {

CCwrapper::CCwrapper() : cc_cut(nullptr), _cutcount(0) {}

CCwrapper::~CCwrapper()
{
  int count = 0;
  lpcut_in *it = cc_cut;
  
  while(it != nullptr){
    ++count;
    cout << "Count is " << count << ", cliq count is : "
	 << it->cliquecount << "\n";
    cc_cut = cc_cut->next;
    CCtsp_free_lpcut_in(it);
    delete(it);
    it = cc_cut;
  }

  cout << "called free/delete " << count << " times\n";
}

}
}
