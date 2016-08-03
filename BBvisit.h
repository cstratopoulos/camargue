#ifndef PSEP_BBVISIT_H
#define PSEP_BBVISIT_H

#include <memory>

#include "BButils.h"
#include "cutcall."
#include "LPprune.h"
#include "LPcore.h"
#include "BBconstraints.h"

namespace PSEP {
  namespace BB {
    class Visitor {
    public:
      Visitor(PSEP::CutControl &_CutControl, PSEP::LPPrune &_LPPrune,
	      PSEP_LP_Core &_LPCore, PSEP::BB::Constraints &_CMgr) :
      CutControl(_CutControl), LPPrune(_LPPrune), LPCore(_LPCore),
	ConstraintMgr(_CMgr) {}

      int previsit(std::unique_ptr<TreeNode> &v);
      int postvisit(std::unique_ptr<TreeNode> &v);

    private:
      int pcut_solve(std::unique_ptr<TreeNode> &v);

      PSEP::CutControl &CutControl;
      PSEP::LPPrune &LPPrune;
      PSEP_LP_Core &LPCore;

      PSEP::BB::Constraints &ConstraintMgr;
    };
  }
}

#endif
