#ifndef CMR_BBVISIT_H
#define CMR_BBVISIT_H

#include <memory>

#include "LPprune.hpp"
#include "LPcore.hpp"
#include "purecut.hpp"
#include "cutcall.hpp"
#include "BButils.hpp"
#include "BBconstraints.hpp"

namespace CMR {
  namespace BB {
    class Visitor {
    public:
      Visitor(CMR::PureCut & _PureCut, CMR::BB::Constraints &_ConsMgr) :
      PureCut(_PureCut), LPPrune(_PureCut.LPPrune), LPCore(_PureCut.LPCore),
	ConstraintMgr(_ConsMgr), RightBranch(ConstraintMgr.RBranch) {}

      int previsit(std::unique_ptr<TreeNode> &v);
      int postvisit(std::unique_ptr<TreeNode> &v);

      //private:
      int handle_augmentation();

      CMR::PureCut &PureCut;
      
      CMR::LP::CutPrune &LPPrune;
      CMR::LP::Core &LPCore;

      CMR::BB::Constraints &ConstraintMgr;
      CMR::BB::RightBranch &RightBranch;
    };
  }
}

#endif
