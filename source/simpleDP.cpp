#include "simpleDP.hpp"
#include "PSEP_util.hpp"

#include <iostream>
#include <memory>

using std::cerr;
using std::unique_ptr;

namespace PSEP {

int Cut<dominoparity>::cutcall()
{
  int rval = 0;

  rval = separate();
  if(rval) goto CLEANUP;

  //rval = add_external

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<dominoparity>::cutcall.\n";
  return rval;
}

int Cut<dominoparity>::separate()
{
  int rval = 0;
  unique_ptr<DPCutGraph> witness;

  rval = candidates.get_light_teeth();
  if(rval) goto CLEANUP;

  rval = candidates.merge_and_sort();
  if(rval) goto CLEANUP;

  candidates.weak_elim();

  try { witness = PSEP::make_unique<DPCutGraph>(candidates); } catch(...){
    PSEP_SET_GOTO(rval, "DPCutGraph allocation failed. ");
  }

  rval = witness->simple_DP_sep(dp_q);
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<dominoparity>::separate.\n";
  return rval;
}

}
