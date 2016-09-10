#include<iostream>
#include<iomanip>

#include<cmath>

#define DO_SAFE_MIR_DBL 1
#include <safemir/src/gen_slvr.hpp>
#include <safemir/src/cplex_slvr.cpp>
#include <safemir/src/ds_slvr.cpp>

#include "safegmi.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

// int Cut<safeGMI>::test(){
//   int frac_basic = -1;
//   cout << "Searching for fractional basic variable...\n";
//   for(int i = 0; i < support_indices.size(); i++){
//     int lp_ind = support_indices[i];
//     if(m_lp_edges[lp_ind] < 1 - LP::EPSILON &&
//        frac_colstat[lp_ind] == CPX_BASIC){
//       frac_basic = lp_ind;
//       cout << "Found fractional basic var: LP entry "
// 	   << m_lp_edges[lp_ind] << ", colstat "
// 	   << frac_colstat[lp_ind] << "\n";
//       break;
//     }
//   }

//   if(frac_basic == -1){
//     cerr << "Somehow couldn't find fractional basic variable, terminating\n";
//     return 1;
//   }

//   cout << "Declaring other parameters for getTableauRow....\n";
//   CUTSsprow_t<double> *tab_row_sp = (CUTSsprow_t<double> *) NULL;
//   int maxnz = PSEPlp_numcols(&m_lp);

//   if(CUTSnewSpRow(&tab_row_sp, maxnz)) return 1;
//   cout << "Declared and allocated sparse tableau row\n";

//   SLVRbasisInfo_t *basis_info = (SLVRbasisInfo_t *) NULL;
//   cout << "Declared null binfo pointer\n";

//   CUTSsystem_t<double> **form_sys = (CUTSsystem_t<double> **) NULL;

//   cout << "Now trying to get tableaurow............";
//   int rval = SLVRgetTableauRow(safemir_LPobj.get(),
// 			       form_sys, &tab_row_sp,
// 			       &basis_info, frac_basic);
//   if(rval)
//     cerr << "Failed\n";
//   else
//     cout << "It worked?\n";

//   cout << "Sense of tableau row: " << char(tab_row_sp->sense) << "\n";

//   CUTSfreeSpRow(&tab_row_sp);
//   CUTSfreeSystem(form_sys);
//   if(basis_info) free(basis_info);
  
//   return 1;
// }

int Cut<safeGMI>::cutcall(){
  cout << "The safe GMI cutcall\n";

  return 1;
}

int Cut<safeGMI>::get_constraint_matrix(){
  int rval = 0;

  rval = PSEPlp_copybase(&m_lp, &frac_colstat[0], &frac_rowstat[0]);
  if(rval){
    fprintf(stderr, "Failed to copy fractional solution basis, ");
    goto CLEANUP;
  }
  
  rval =  PSEPlp_no_opt(&m_lp);
  if(rval){
    fprintf(stderr, "Failed to pivot back to fractional solution, ");
    goto CLEANUP;
  }

  try { safe_mir_data.reset(new SafeMIRGroup(m_lp)); }
  catch(const std::bad_alloc &e){
    fprintf(stderr, "Out of memory for new safe_mir_data, ");
    rval = 1; goto CLEANUP;
  }

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<safeGMI>::get_constraint_matrix\n";
  return rval;
}
