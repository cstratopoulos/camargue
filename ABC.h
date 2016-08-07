#ifndef PSEP_ABC_H
#define PSEP_ABC_H

#include<memory>
#include<stack>

#include "BButils.h"
#include "BBconstraints.h"
#include "BBvisit.h"
#include "datagroups.h"
#include "cutcall.h"

namespace PSEP {
  class ABC {
  public:
  ABC(PSEP_GraphGroup &GraphGroup, PSEP_BestGroup &BestGroup,
      PSEP_LPGroup &LPGroup, PSEP_SupportGroup &SupportGroup,
      std::vector<double> &lower_bounds,
      PSEP::PureCut &_PureCut):
    EdgeStats(lower_bounds),
      ConstraintMgr(GraphGroup, BestGroup, LPGroup, SupportGroup,
		    RightBranch, EdgeStats, _PureCut.LPPrune,
		    _PureCut.LPcore),
      Visitor(_PureCut, ConstraintMgr){}
  
    int solve();
    
  private:
    std::stack<std::unique_ptr<PSEP::BB::TreeNode> > BBTree;
    PSEP::BB::EdgeStatuses EdgeStats;
    PSEP::BB::RightBranch RightBranch;
    PSEP::BB::Constraints ConstraintMgr;
    PSEP::BB::Visitor Visitor;
  };
}

#endif
