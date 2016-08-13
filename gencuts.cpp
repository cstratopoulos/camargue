#include<iostream>

#include<cmath>

#include "gencuts.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

int Cut<general>::separate(const double piv_val){
  int rval = 0;
  int num_frac, num_mir, num_disj, num_total = 0;
  int numcols = PSEPlp_numcols(&m_lp), numrows = PSEPlp_numrows(&m_lp);
  double objval;
  mip_cut_candidates generated_cuts(numcols, numrows);

  rval = init_mip(piv_val, generated_cuts);
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

  rval = check_cuts(generated_cuts);
  if(rval) goto CLEANUP;


  rval = revert_lp();
  if(rval) goto CLEANUP;

  if(num_total == 0)
    rval = 2;

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<general>::separate\n";
  return rval;
}

int Cut<general>::init_mip(const double piv_val,
			   mip_cut_candidates &callback_args){
  int numrows = PSEPlp_numrows(&m_lp);
  double cutfactor = (numrows + (double) max_cuts) / numrows;

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

  rval = CPXsetbranchcallbackfunc(m_lp.cplex_env, &branchcallback,
				  &callback_args);
  if(rval) { cerr << "Couldn't set callback func, "; goto CLEANUP; }

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

  rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_Threads, 1);
  if(rval) { cerr << "Couldn't clamp threads, "; goto CLEANUP; }

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

  rval = CPXsetintparam(m_lp.cplex_env, CPXPARAM_Threads, 0);
  if(rval) { cerr << "Couldn't revert threads, "; goto CLEANUP; }

  rval = CPXsetbranchcallbackfunc(m_lp.cplex_env, NULL, NULL);
  if(rval) { cerr << "Couldn't get rid of branch callback, "; goto CLEANUP; }

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

//callback code adapted from
// https://www.ibm.com/developerworks/community/forums/html/
// topic?id=77777777-0000-0000-0000-000014468982

int Cut<general>::branchcallback (CPXCENVptr xenv, void *cbdata, int wherefrom,
			   void *cbhandle, int brtype, int brset, int nodecnt,
			   int bdcnt, const int *nodebeg, const int *xindex,
			       const char *lu, const double *bd,
			       const double *nodeest, 
			       int *useraction_p){
  (void) brtype; (void) brset; (void) nodecnt; (void) bdcnt; (void) nodebeg;
  (void) xindex; (void) lu; (void) bd; (void) useraction_p;


  generated_cut *const arg = (generated_cut *) cbhandle;
  CPXLPptr nodelp;
  int rval = 0, rows, cols;

  rval = CPXgetcallbacknodelp(xenv, cbdata, wherefrom, &nodelp);
  if(rval) { fprintf(stderr, "CPXgetcallbacknodelp failed, "); goto CLEANUP;}
  rows = CPXgetnumrows(xenv, nodelp);
  cols = CPXgetnumcols(xenv, nodelp);


    for(int i = arg->next_cut; i < rows; i++, j++){
      int rmatbeg, surplus, nzcount;
      char rhs;

      rval = CPXgetrows(xenv, nodelp, &nzcount, &rmatbeg,
			&(arg->index_buffer[0]),
			&(arg->coefficient_buffer[0]), cols,
			&surplus, i, i);
      if(rval) { fprintf(stderr, "CPXgetrows failed, "); goto CLEANUP; }

      rval = CPXgetsense(xenv, nodelp, &(arg->senses)[j], i, i);
      if(rval) { fprintf(stderr, "CPXgetsense failed, "); goto CLEANUP; }

      rval = CPXgetrhs(xenv, nodelp, &(arg->rhs_array)[j], i, i);
      if(rval) { fprintf(stderr, "CPXgetrhs failed, "); goto CLEANUP; }
    }
    

  if(rows > arg->next_cut)
    arg->next_cut = rows;

 CLEANUP:
  if(rval)
    fprintf(stderr, "problem in Cut<general>::branch_callback, rval %d\n",
	    rval);
  return rval;
}
