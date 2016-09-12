#include<algorithm>
#include<iostream>
#include<iomanip>

#include<cmath>

#define DO_SAFE_MIR_DBL 1
#define SAFE_MIR_DEBUG_LEVEL DBG_LEVEL_HIGH
#define CUTSslackSign( row ) ( row->sense == 'L' ? 1 : (row->sense == 'E' ? 0 : -1 ) )
#include<safemir/src/gen_slvr.hpp>
#include<safemir/src/util_cuts.hpp>
#include<safemir/src/gen_mir.cpp>
#include<safemir/src/cplex_slvr.cpp>
#include<safemir/src/ds_slvr.cpp>
#include<safemir/src/ds_cuts.cpp>
#include<safemir/src/safe_mir_dbl.cpp>

#include "safegmi.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

int Cut<safeGMI>::cutcall(){
  int rval = 0;
  
  cout << "The safe GMI cutcall\n";
  rval = init_constraint_info();
  if(rval) goto CLEANUP;

  rval = get_cuts();
  if(rval) goto CLEANUP;
  else cout << "Successful call\n";

  
 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<safeGMI>::cutcall\n";
  safe_mir_data.reset(NULL);
  return 1;
}

int Cut<safeGMI>::init_constraint_info(){
  int rval = 0;
  int numcols = m_lp_edges.size();
  int numrows = PSEPlp_numrows(&m_lp);

  cout << "LP currently has " << numcols << " cols, " << numrows << " rows, "
       << (numcols + numrows) << " total\n";

  rval = PSEPlp_copystart(&m_lp, &frac_colstat[0], &frac_rowstat[0],
			  &m_lp_edges[0], NULL, NULL, NULL);
  if(rval) GOTO_CLEANUP("Failed to copy frac solution, ");

  rval = PSEPlp_no_opt(&m_lp);
  if(rval) GOTO_CLEANUP("Failed to factor basis, ");
  
  cout << "Copied frac basis and pivoted back\n";
  
  try { safe_mir_data.reset(new SafeMIRGroup(m_lp)); }
  catch(const std::bad_alloc &e){
    rval = 1; GOTO_CLEANUP("Out of memory for safe_mir_data reset, ");
  }
  cout << "Reset data unique_ptr\n";

  rval = SLVRformulationRows(&(safe_mir_data->lp_obj),
			     &(safe_mir_data->constraint_matrix));
  if(rval) GOTO_CLEANUP("SLVRformulationRows failed, ");

  cout << "Got formulation rows\n";

  {//these braces give form_sys pointer local scope; it's a naming convenience
    CUTSsystem_t<double> *form_sys = safe_mir_data->constraint_matrix;

    for(int i = 0; i < form_sys->sys_rows; i++){
      rval = CUTSaddSlackVariable(form_sys->rows[i], numcols + i);
      if(rval) GOTO_CLEANUP("CUTSaddSlackVariable failed, ");
    }
  }
  cout << "Added slack variables\n";

  rval = SLVRgetBasisInfo(&(safe_mir_data->lp_obj),
			  &(safe_mir_data->basis_info));
  if(rval) GOTO_CLEANUP("SLVRgetBasisInfo failed, ");
  cout << "Got basis info\n";

  rval = SLVRgetVarInfo(&(safe_mir_data->lp_obj), true,
			&(safe_mir_data->var_info));
  if(rval) GOTO_CLEANUP("SLVRgetVarInfo failed, ");
  cout << "Got var info, with slacks\n";

  safe_mir_data->full_x = SLVRgetFullX(&(safe_mir_data->lp_obj),
				       safe_mir_data->constraint_matrix,
				       &m_lp_edges[0]);
  if(!safe_mir_data->full_x) GOTO_CLEANUP("SLVRgetFullX failed, ");
  cout << "Got LP solution plus slacks (full x)\n";

  rval = CUTSnewFlips(&(safe_mir_data->flips), numcols + numrows);
  if(rval) GOTO_CLEANUP("CUTSnewFlips failed, ");

  rval = MIRwhoToFlip(safe_mir_data->var_info,
		      safe_mir_data->full_x,
		      safe_mir_data->flips, numcols + numrows);
  if(rval) GOTO_CLEANUP("MIRwhoToFlip failed, ");
  cout << "Initialized flips and retrieved flip info, flips->nvars:"
       << safe_mir_data->flips->nvars << "\n";
  

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<safeGMI>::get_constraint_matrix\n";
  return rval;
}

int Cut<safeGMI>::get_cuts(){
  int rval = 0;
  int numcols = m_lp_edges.size();
  int numrows = PSEPlp_numrows(&m_lp);

  vector<int> frac_basic_vars;
  for(int i = 0; i < support_indices.size(); i++){
    int index = support_indices[i];
    if((m_lp_edges[index] < 1 - LP::EPSILON) &&
       frac_colstat[index] == CPX_BASIC){
      try{ frac_basic_vars.push_back(index); }
      catch(const std::bad_alloc &){
	rval = 1; GOTO_CLEANUP("Out of memory for frac_basic_vars, ");
      }
    } 
  }

  sort(frac_basic_vars.begin(), frac_basic_vars.end(),
       [](const double &x, const double &y) -> bool {
	 return fabs(x - 0.5) < fabs(y - 0.5);
       });

  cout << "Created and sorted list of fractional basic vars\n";

  rval = (CUTSnewSpRow(&(safe_mir_data->tab_row_sparse),
		       numcols + numrows) ||
	  CUTSnewSpRow(&(safe_mir_data->current_cut_sparse),
		       numcols + numrows) ||
	  CUTSnewSpRow(&(safe_mir_data->best_cut_sparse),
		       numcols + numrows));
  if(rval) GOTO_CLEANUP("Out of memory for tab/cut rows, ");
  cout << "Allocated empty sparse cut/tab rows\n";
  
  for(int i = 0; i < frac_basic_vars.size(); i++){
    int current_var = frac_basic_vars[i];
    double modified_rhs = 0;
    rval = SLVRgetTableauRow(&(safe_mir_data->lp_obj),
			     &(safe_mir_data->constraint_matrix),
			     &(safe_mir_data->tab_row_sparse),
			     &(safe_mir_data->basis_info), current_var);
    if(rval) GOTO_CLEANUP("SLVRgetTableauRow failed, ");
    cout << "Tab row " << i << " has sense "
	 << safe_mir_data->tab_row_sparse->sense << ", "
	 << safe_mir_data->tab_row_sparse->nz << " nonzeros, rhs: "
	 << safe_mir_data->tab_row_sparse->rhs << ", with slack: "
	 << safe_mir_data->tab_row_sparse->with_slack << "\n";

    // problem here: supposedly the tableau row contains a higher index
    // than numcols. need to check if it contains slacks somehow?
    rval = MIRsafeComputeModRhs_dbl(safe_mir_data->tab_row_sparse,
    				    safe_mir_data->var_info,
    				    &modified_rhs,
    				    safe_mir_data->flips);
    if(rval) GOTO_CLEANUP("MIRsafeComputeModRhs_dbl failed, ");
    cout << "Modified rhs: " << modified_rhs << "\n";
				    
  }
  cout << "Got all tab rows and did nothing with them\n";

 CLEANUP:
  if(rval == 1)
    cerr << "problem in Cut<safeGMI>::get_cuts\n";
  return rval;
}
