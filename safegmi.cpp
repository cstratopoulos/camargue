#include<iostream>
#include<iomanip>

#include<cmath>

#define DO_SAFE_MIR_DBL 1
#define CUTSslackSign( row ) ( row->sense == 'L' ? 1 : (row->sense == 'E' ? 0 : -1 ) )
#include <safemir/src/gen_slvr.hpp>
#include <safemir/src/util_cuts.hpp>
#include <safemir/src/cplex_slvr.cpp>
#include <safemir/src/ds_slvr.cpp>

#include "safegmi.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

int Cut<safeGMI>::cutcall(){
  int rval = 0;
  
  cout << "The safe GMI cutcall\n";
  rval = get_constraint_matrix();
  if(rval) goto CLEANUP;
  else cout << "Successful call\n";

  
 CLEANUP:
  safe_mir_data.reset(NULL);
  return 1;
}

int Cut<safeGMI>::get_constraint_matrix(){
  int rval = 0;

  rval = PSEPlp_copybase(&m_lp, &frac_colstat[0], &frac_rowstat[0]);
  if(rval){
    cerr << "Failed to copy fractional solution basis, ";
    goto CLEANUP;
  }
  
  rval =  PSEPlp_no_opt(&m_lp);
  if(rval){
    cerr << "Failed to pivot back to fractional solution, ";
    goto CLEANUP;
  }
  cout << "Copied frac basis and pivoted back\n";

  try { safe_mir_data.reset(new SafeMIRGroup(m_lp)); }
  catch(const std::bad_alloc &e){
    cerr << "Out of memory for new safe_mir_data, ";
    rval = 1; goto CLEANUP;
  }
  cout << "Reset data unique_ptr\n";

  rval = SLVRformulationRows(&(safe_mir_data->lp_obj),
			     &(safe_mir_data->constraint_matrix));
  if(rval){
    cerr << "SLVRformulationRows failed, "; goto CLEANUP;
  }

  cout << "Got formulation rows\n";

  {
    CUTSsystem_t<double> *form_sys = safe_mir_data->constraint_matrix;
    int ncols = PSEPlp_numcols(&m_lp);

    for(int i = 0; i < form_sys->sys_rows; i++){
      rval = CUTSaddSlackVariable(form_sys->rows[i], ncols + i);
      if(rval){
	cerr << "CUTSaddSlackVariable failed, "; goto CLEANUP;
      }
    }
  }

  cout << "Added slack variables\n";

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<safeGMI>::get_constraint_matrix\n";
  return rval;
}
