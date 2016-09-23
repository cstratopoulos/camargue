#include<iostream>
#include<vector>

#include "BBvisit.hpp"
#include "PSEP_util.hpp"
#include "pivplan.hpp"

using namespace std;
using namespace PSEP;
using namespace PSEP::BB;

int Visitor::previsit(unique_ptr<TreeNode> &v){
  if(v->type() == NodeType::ROOT){
    v->node_stat = NodeStat::VISITED;
    return 0;
  }
  
  int rval = 0;
  LP::PivType piv_status;
  PivotPlan plan(ConstraintMgr.ncount, PivPresets::BRANCH);

  rval = ConstraintMgr.enforce(v);
  if(rval) goto CLEANUP;

  if(v->type() == NodeType::RIGHT){
    rval = ConstraintMgr.prune();
    if(rval) goto CLEANUP;
  }

  while(true){
    rval = PureCut.solve(plan, piv_status);

    if(piv_status == LP::PivType::Tour){
      rval = handle_augmentation();
      if(rval) goto CLEANUP;
      continue;
    }

    break;
  }

  switch(piv_status){
  case LP::PivType::FathomedTour:
    v->node_stat = NodeStat::FATHOMED;
    break;
  // case PivType::SUBTOUR:
  //   cerr << "Pivoted to subtour with no segment cut found!\n";
  //   rval = 1;
  //   goto CLEANUP;
  default:
    v->node_stat = NodeStat::VISITED;
    break;
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in Visitor::previsit\n";
  return rval;
}

int Visitor::postvisit(unique_ptr<TreeNode> &v){
  int rval = ConstraintMgr.unenforce(v);

  if(rval)
    cerr << "Problem in Visitor::postvisit\n";
  return rval;
}

int Visitor::handle_augmentation(){
  int rval = 0;

  if(RightBranch.active()){
    rval = ConstraintMgr.update_right_rows();
    if(rval) goto CLEANUP;
  }

  rval = LPCore.update_best_tour();
  if(rval) goto CLEANUP;

  rval = ConstraintMgr.prune();
  if(rval) goto CLEANUP;


 CLEANUP:
  if(rval)
    cerr << "Problem in Visitor::handle_augmentation\n";
  return rval;
}
