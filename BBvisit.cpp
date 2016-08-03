#include<iostream>

#include "BBvisit.h"

using namespace std;
using namespace PSEP::BB;

int Visitor::previsit(unique_ptr<TreeNode> &v){
  int rval = ConstraintMgr.enforce(v);
  if(rval) goto CLEANUP;

  rval = pcut_solve(v);

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
