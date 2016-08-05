#include "ABC.h"

using namespace std;
using namespace PSEP;
using namespace PSEP::BB;

int ABC::solve(){
  int rval = 0;
  BBTree.push(unique_ptr<TreeNode>(new TreeNode));

  while(!BBTree.empty()){
    if(BBTree.top()->status() == NodeStat::UNVISITED){
      rval = Visitor.previsit(BBTree.top());
      if(rval) goto CLEANUP;

      if(BBTree.top()->status() != NodeStat::FATHOMED){
	int newedge = ConstraintMgr.compute_branch_edge();
	rval = (newedge == -1);
	if(rval){
	  cerr << "No new branch edge found\n"; goto CLEANUP;
	}

	BBTree.push(unique_ptr<TreeNode>(new TreeNode(NodeType::RIGHT,
							   newedge)));
	BBTree.push(unique_ptr<TreeNode>(new TreeNode(NodeType::LEFT,
							   newedge)));
	continue;
      }
    }

    rval = Visitor.postvisit(BBTree.top());
    if(rval) goto CLEANUP;
    BBTree.pop();
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in ABC::solve\n";
  return rval;
}
