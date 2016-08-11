#include<iostream>

#include<cmath>

#include "gencuts.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

int Cut<general>::separate(const double piv_val){
  int rval = 0;
  int num_frac, num_mir, num_disj, num_total = 0;
  double objval;

  rval = init_mip(piv_val);
  if(rval) goto CLEANUP;

  rval = make_all_binary();
  if(rval) goto CLEANUP;

  rval = PSEPmip_opt(&m_lp);
  if(rval) goto CLEANUP;

  rval = PSEPmip_getbestobjval(&m_lp, &objval);
  if(rval) goto CLEANUP;

  cout << "    Mip opt has objval: " << objval << "\n";

  rval = num_added(num_frac, num_disj, num_mir);
  num_total = num_frac + num_disj + num_mir;
  if(rval) goto CLEANUP;
  cout << "\n    Added " << num_frac << " gomory fractional cuts, "
       << num_disj << " disjunctive cuts, "
       << num_mir << " MIR cuts\n\n";


  rval = revert_lp();
  if(rval) goto CLEANUP;

  if(num_total == 0)
    rval = 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<general>::separate\n";
  return rval;
}

int Cut<general>::init_mip(const double piv_val){
  int numrows = PSEPlp_numrows(&m_lp);
  double cutfactor = (numrows + 10.0) / numrows;

  deletion_row = numrows;
  int numcols = PSEPlp_numcols(&m_lp);
  int rmatbeg = 0;
  char sense = 'G';
  double rhs = piv_val;
  
  int rval = 0;

  vector<int> lp_indices(numcols);
  vector<double> obj_fun(numcols);
  for(int i = 0; i < numcols; i++) lp_indices[i] = i;

  rval = PSEPlp_make_mip(&m_lp);
  if(rval) goto CLEANUP;

  rval = PSEPlp_getobj(&m_lp, &obj_fun[0], numcols);
  if(rval) goto CLEANUP;

  cout << "Making objective function >= " << rhs << "....";
  rval = PSEPlp_addrows(&m_lp, 1, numcols, &rhs, &sense, &rmatbeg,
			&lp_indices[0], &obj_fun[0]);
  if(rval) goto CLEANUP;
  cout << "Done\n";

  if(gencuts.gomory_frac){
    rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_MIP_Cuts_Gomory, 2);
    if(rval){
      cerr << "Couldn't enable Gomory cuts, ";
      goto CLEANUP;
    }
  }

  if(gencuts.disjunctive){
    rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_MIP_Cuts_Disjunctive, 3);
    if(rval){
      cerr << "Couldn't enable disjunctive cuts, ";
      goto CLEANUP;
    }
  }

  if(gencuts.rounding){
    rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_MIP_Cuts_MIRCut, 1);
    if(rval){
      cerr << "Couldn't enable MIR cuts, ";
      goto CLEANUP;
    }
  }

  rval = CPXsetdblparam(m_lp.cplex_env, CPXPARAM_MIP_Limits_CutsFactor,
  			cutfactor);
  if(rval){
    cerr << "Couldn't set cut row factor, ";
    goto CLEANUP;
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in Cut<general>::init_mip\n";
  return rval;
}

int Cut<general>::revert_lp(){
  int rval = PSEPlp_make_lp(&m_lp);
  if(rval){
    cerr << "Couldn't change back to LP, ";
    goto CLEANUP;
  }

  cout << "Removing objective function bound....";
  rval = PSEPlp_delrows(&m_lp, deletion_row, deletion_row);
  if(rval) goto CLEANUP;
  cout << "Done.\n";

  rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_MIP_Cuts_Gomory, -1);
  if(rval){
    cerr << "Couldn't disable Gomory cuts, ";
    goto CLEANUP;
  }

  rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_MIP_Cuts_Disjunctive, -1);
  if(rval){
    cerr << "Couldn't disable disjunctive cuts, ";
    goto CLEANUP;
  }
  
  rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_MIP_Cuts_MIRCut, -1);
  if(rval){
    cerr << "Couldn't disable MIR cuts, ";
    goto CLEANUP;
  }

 CLEANUP:
  if(rval)
    cerr << "Problem in Cut<general>::revert_lp\n";
  return rval;
}

int Cut<general>::make_all_binary(){
  int numcols = PSEPlp_numcols(&m_lp);
  vector<char> vartype(numcols, CPX_BINARY);
  vector<int> indices(numcols);
  
  for(int i = 0; i < numcols; i++)
    indices[i] = i;

  int rval = PSEPmip_change_vartype(&m_lp, numcols, &indices[0], &vartype[0]);
  if(rval)
    cerr << "Problem in Cut<general>::make_all_binary\n";
  return rval;
}

int Cut<general>::make_binary(const int edge){
  char vartype = CPX_BINARY;

  int rval = PSEPmip_change_vartype(&m_lp, 1, &edge, &vartype);
  if(rval)
    cerr << "Problem in Cut<general>::make_binary\n";
  return rval;
}

int Cut<general>::num_added(int &frac, int &disj, int &mir){
  int rval = PSEPmip_getnumcuts(&m_lp, CPX_CUT_FRAC, &frac);
  if(rval) goto CLEANUP;

  rval = PSEPmip_getnumcuts(&m_lp, CPX_CUT_DISJ, &disj);
  if(rval) goto CLEANUP;
  
  rval = PSEPmip_getnumcuts(&m_lp, CPX_CUT_MIR, &mir);
  if(rval) goto CLEANUP;

 CLEANUP:
  if(rval)
    cerr << "Problem in Cut<general>::num_added\n";
  return rval;
}