#ifndef PSEP_BBVISIT_H
#define PSEP_BBVISIT_H

#include <memory>

#include "LPprune.h"
#include "LPcore.h"
#include "purecut.h"
#include "cutcall.h"
#include "BButils.h"
#include "BBconstraints.h"

namespace PSEP {
  namespace BB {
    class Visitor {
    public:
      Visitor(PSEP::PureCut & _PureCut, PSEP::BB::Constraints &_ConsMgr) :
      PureCut(_PureCut), LPPrune(_PureCut.LPPrune), LPCore(_PureCut.LPCore),
	ConstraintMgr(_ConsMgr), RightBranch(ConstraintMgr.RBranch) {}

      int previsit(std::unique_ptr<TreeNode> &v);
      int postvisit(std::unique_ptr<TreeNode> &v);

      //private:
      int handle_augmentation();

      PSEP::PureCut &PureCut;
      
      PSEP::LP::CutPrune &LPPrune;
      PSEP::LP::Core &LPCore;

      PSEP::BB::Constraints &ConstraintMgr;
      PSEP::BB::RightBranch &RightBranch;
    };
  }
}

#endif
