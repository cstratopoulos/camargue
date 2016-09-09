#include<iostream>
#include<iomanip>
#include<memory>

#include<cmath>

#include<safemir/src/gen_slvr.hpp>
#include<safemir/src/cplex_slvr.hpp>


#include "gen_slvr.hpp"
#include "cplex_slvr.hpp"


#include "safegmi.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

int Cut<safeGMI>::test(){
  cout << "Testing safeGMI procedures\n";
  double obj;
  PSEPlp_objval(&m_lp, &obj);
  cout << "Objval supposedly: " << obj << "\n";

  cout << "Trying to pivot back to frac solution....\n";
  PSEPlp_copybase(&m_lp, &frac_colstat[0], &frac_rowstat[0]);
  PSEPlp_no_opt(&m_lp);
  PSEPlp_objval(&m_lp, &obj);
  cout << "Objval now: " << obj << "\n";

  cout << "Creating vartype vector...\n";
  std::vector<char> vartype(PSEPlp_numcols(&m_lp), 'B');

  cout << "Attempting to create safemir cplex object....";
  unique_ptr<SLVRcplex_t> safemir_LPobj(new SLVRcplex_t);
  safemir_LPobj->env = m_lp.cplex_env;
  safemir_LPobj->prob = m_lp.cplex_lp;
  safemir_LPobj->ctype = &vartype[0];
  cout << "Done?\n";

  int frac_basic = -1;
  cout << "Searching for fractional basic variable...\n";
  for(int i = 0; i < support_indices.size(); i++){
    int lp_ind = support_indices[i];
    if(m_lp_edges[lp_ind] < 1 - LP::EPSILON &&
       frac_colstat[lp_ind] == CPX_BASIC){
      frac_basic = lp_ind;
      cout << "Found fractional basic var: LP entry "
	   << m_lp_edges[lp_ind] << ", colstat "
	   << frac_colstat[lp_ind] << "\n";
      break;
    }
  }

  if(frac_basic == -1){
    cerr << "Somehow couldn't find fractional basic variable, terminating\n";
    return 1;
  }

  cout << "Declaring other parameters for getTableauRow....\n";
  CUTSsprow_t<double> *tab_row_sp = (CUTSsprow_t<double> *) NULL;
  int maxnz = PSEPlp_numcols(&m_lp);

  if(CUTSnewSpRow(&tab_row_sp, maxnz)) return 1;
  cout << "Declared and allocated sparse tableau row\n";

  SLVRbasisInfo_t *basis_info = (SLVRbasisInfo_t *) NULL;
  cout << "Declared null binfo pointer\n";

  CUTSsystem_t<double> **form_sys = (CUTSsystem_t<double> **) NULL;

  cout << "Now trying to get tableaurow............";
  int rval = SLVRgetTableauRow(safemir_LPobj.get(),
			       form_sys, &tab_row_sp,
			       &basis_info, frac_basic);
  if(rval)
    cerr << "Failed\n";
  else
    cout << "It worked?\n";

  CUTSfreeSpRow(&tab_row_sp);
  if(basis_info) free(basis_info);
  
  return 1;
}
