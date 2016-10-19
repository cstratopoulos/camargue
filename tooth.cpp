#include "tooth.hpp"

#include <iostream>

using namespace PSEP;
using namespace std;

inline bool SimpleTooth::sandwich() const
{
  return
    (body_start <= body_end) ?
    (body_start <= root && root <= body_end) : //[___<----*--->__]
    (body_start <= root || root <= body_end); // [-->__<--*--] OR [-*->__<--]
}

inline int CandidateTeeth::body_size(const SimpleTooth &T)
{
  return
    (T.body_start <= T.body_end) ?
    (T.body_end - T.body_start + !T.sandwich()) :
    ((G_s.node_count - T.body_start) + T.body_end + !T.sandwich());    
}

void CandidateTeeth::complement(SimpleTooth &T)
{
  int ncount = G_s.node_count;
  int body_start = T.body_start, body_end = T.body_end, root = T.root;
  int c_start, c_end;

  if((body_start == ((root + 1) % ncount)) ||
     (root == ((body_end + 1) % ncount))){
    if(body_start == ((root + 1) % ncount)){
      c_start = (body_end + 1) % ncount;
      c_end = (root + ncount - 1) % ncount;
    } else {
      c_start = (root + 1) % ncount;
      c_end = (body_end + ncount - 1) % ncount;
    }
  } else {
    c_start = (body_end + 1) % ncount;
    c_end = (body_start + ncount - 1) % ncount;
  }

  T.body_start = c_start;
  T.body_end = c_end;
}

int CandidateTeeth::body_subset(const SimpleTooth &T, const SimpleTooth &R,
				bool &result)
{
  if(T.root != R.root){
    cerr << "Cannot currently test body subset with different roots\n";
    return 1;
  }

  if(R.body_start <= R.body_end){ //if R is contiguous
    if(T.body_start > T.body_end){//R cannot contain a sandwich tooth T
      result = false;
    } else {
      result = (R.body_start <= T.body_start && T.body_end <= R.body_end);
    }
  } else { //if R wraps
    result = ( (R.body_start <= T.body_start &&
		(R.body_start <= T.body_end || T.body_end <= R.body_end)) ||
	       (T.body_start <= R.body_end && T.body_end <= R.body_end));
  }

  return 0;
}
