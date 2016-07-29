#ifndef PSEP_CUTS_H
#define PSEP_CUTS_H

#include <vector>
#include <memory>

#include "datagroups.h"
#include "Graph.h"

namespace PSEP {
  template<class cut_t>
    class Cuts {
  public:
  Cuts(PSEP_GraphGroup &_GraphGroup, PSEP_BestGroup &_BestGroup,
       PSEP_SupportGroup &_SupportGroup, PSEP_LPGroup &_LPGroup):
    GraphGroup(_GraphGroup), BestGroup(_BestGroup), SupportGroup(_SupportGroup),
      LPGroup(_LPGroup){}
    
    int separate(){return 1;}
    
    int get_coefficients(int *deltacount){*deltacount = -1; return 1;}
 
    int add_cut(const int deltacount){return 1;}
    int add_cut(const int deltacount, const int cutedge){return 1;};

    int cut_call(){return 1;}

  private:
    PSEP_GraphGroup &GraphGroup;
    PSEP_BestGroup &BestGroup;
    PSEP_SupportGroup &SupportGroup;
    PSEP_LPGroup &LPGroup;

    std::unique_ptr<cut_t> best;
  };
}

#endif
