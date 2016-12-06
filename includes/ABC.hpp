#ifndef CMR_ABC_H
#define CMR_ABC_H

#include<memory>
#include<stack>

#include "datagroups.hpp"
#include "cutcall.hpp"
#include "BButils.hpp"
#include "BBconstraints.hpp"
#include "BBvisit.hpp"

namespace CMR {
  class ABC {
  public:
  ABC(CMR::BB::BranchPlan Strategy,
      Data::GraphGroup &GraphGroup, Data::BestGroup &BestGroup,
      Data::LPGroup &LPGroup, Data::SupportGroup &SupportGroup,
      std::vector<double> &lower_bounds,
      CMR::PureCut &_PureCut):
    EdgeStats(lower_bounds),
      ConstraintMgr(Strategy, GraphGroup, BestGroup, LPGroup, SupportGroup,
		    RightBranch, EdgeStats, _PureCut.LPPrune,
		    _PureCut.LPCore),
      Visitor(_PureCut, ConstraintMgr){}
  
    int solve();
    
  private:
    std::stack<std::unique_ptr<CMR::BB::TreeNode> > BBTree;
    CMR::BB::EdgeStatuses EdgeStats;
    CMR::BB::RightBranch RightBranch;
    CMR::BB::Constraints ConstraintMgr;
    CMR::BB::Visitor Visitor;
  };
}

#endif
