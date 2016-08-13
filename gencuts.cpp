#include<iostream>

#include<cmath>

#include "gencuts.h"
#include "mip.h"

using namespace std;
using namespace PSEP;

int Cut<general>::separate(const double piv_val){
  int rval = 0;
  int numcols = PSEPlp_numcols(&m_lp), numrows = PSEPlp_numrows(&m_lp);
  int rmatbeg = 0, num_nonzeros = 0;
  double objval;
  generated_cut candidate(numcols, numrows + 1, best_tour_edges, m_lp_edges);

  double lhs_lp = 0, lhs_best = 0;

  rval = init_mip(piv_val, candidate);
  if(rval) goto CLEANUP;

  rval = make_all_binary();
  if(rval) goto CLEANUP;
  
  rval = PSEPmip_opt(&m_lp);
  if(rval) goto CLEANUP;

  rval = PSEPmip_getbestobjval(&m_lp, &objval);
  if(rval) goto CLEANUP;

  rval = revert_lp();
  if(rval) goto CLEANUP;

  num_nonzeros = candidate.best_nonzeros;

  if(num_nonzeros == 0){
    rval = 2; goto CLEANUP;
  }

  cout << "Printing best stuff...is tight: "
       << candidate.is_best_exact << "\n";
  cout << "num nz: " << num_nonzeros << "\n";
  // cout << "Indices (size " << candidate.best_indices.size() << ")\n";
  // for(int i = 0; i < candidate.best_indices.size(); i++){
  //   if(i == num_nonzeros) cout << "**cutoff**\n";
  //   cout << candidate.best_indices[i] << endl;   
  // }

  // cout << "Coeffs:\n";
  // for(int i = 0; i < candidate.best_coeffs.size(); i++){
  //   if(i == num_nonzeros) cout << "**cutoff**\n";
  //   cout << candidate.best_coeffs[i] << endl;
  // }

  cout << "rhs: " << candidate.best_rhs << "\n";
  cout << "viol: " << candidate.best_viol << "\n";

  cout << "Double checking by struct...";
  for(int i = 0; i < num_nonzeros; i++){
    int index = candidate.best_indices[i];
    double coeff = candidate.best_coeffs[i];

    lhs_best += best_tour_edges[index] * coeff;
    lhs_lp += m_lp_edges[index] * coeff;
  }

  cout << "Is tight: " << (lhs_best  == candidate.best_rhs) << "\n";
  cout << "lp lhs: " << lhs_lp << "\n";
  cout << "viol: " << fabs(lhs_lp - candidate.best_rhs) << "\n";

  return 1;
    

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<general>::separate\n";
  return rval;
}

int Cut<general>::init_mip(const double piv_val,
			   generated_cut &callback_arg){
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

  cout << "Setting branch callback func...";
  rval = CPXsetsolvecallbackfunc(m_lp.cplex_env, &solvecallback,
				  &callback_arg);
  if(rval) { cerr << "Couldn't set callback func, "; goto CLEANUP; }
  cout << "Done\n";

  rval = PSEPlp_getobj(&m_lp, &obj_fun[0], numcols);
  if(rval) goto CLEANUP;

  rval = PSEPlp_addrows(&m_lp, 1, numcols, &rhs, &sense, &rmatbeg,
			&lp_indices[0], &obj_fun[0]);
  if(rval) goto CLEANUP;

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

  rval = CPXsetsolvecallbackfunc(m_lp.cplex_env, NULL, NULL);
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

int CPXPUBLIC Cut<general>::solvecallback(CPXCENVptr env, void *cbdata,
					  int wherefrom, void *cbhandle,
					  int *useraction_p){
  int rval = 0, numrows, numcols;
  CPXLPptr nodelp;
  generated_cut *arg = (generated_cut *)cbhandle;

  *useraction_p = CPX_CALLBACK_DEFAULT;
  printf("   CALLBACK INVOCATION -- ");

  rval = CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp);
  if(rval) { fprintf(stderr, "CPXgetcallbacknodelp"); goto CLEANUP; }

  numrows = CPXgetnumrows(env, nodelp);
  numcols = CPXgetnumcols(env, nodelp);
  printf(" %d rows in LP, %d initial\n", numrows, arg->initial_numrows);

  if(numrows > arg->initial_numrows){
    int rmatbeg, surplus, num_nonzero;
    char sense;
    double rhs;
    for(int i = arg->initial_numrows; i < numrows; i++){
      rval = CPXgetrows(env, nodelp, &num_nonzero, &rmatbeg,
			&(arg->index_buffer[0]),
			&(arg->coefficient_buffer[0]), numcols, &surplus, i, i);
      if(rval) { fprintf(stderr, "CPXgetrows"); goto CLEANUP; }

      rval = CPXgetsense(env, nodelp, &sense, i, i);
      if(rval) { fprintf(stderr, "CPXgetsense"); goto CLEANUP; }

      rval = CPXgetrhs(env, nodelp, &rhs, i, i);
      if(rval) { fprintf(stderr, "CPXgetrhs"); goto CLEANUP; }

      int ind;
      double coeff;
      double lhs_lp = 0, lhs_best = 0;
      printf("ANALYZING ROW NUMBER %d\n", i);
      
      for(int j = 0; j < num_nonzero; j++){		
	ind = arg->index_buffer[j];
	coeff = arg->coefficient_buffer[j];
	lhs_lp += arg->m_lp_edges[ind] * coeff;
	lhs_best += arg->best_tour_edges[ind] * coeff;
      }
      bool lp_violated;
      bool feasible;
      double viol;
      bool exact_cut = lhs_best == rhs;
      bool epsilon_feasible;

      switch(sense){
      case 'L':
	lp_violated = lhs_lp > rhs;
	viol = lhs_lp - rhs;
	feasible = lhs_best <= rhs;
	break;
      case 'G':
	lp_violated = lhs_lp < rhs;
	viol = rhs - lhs_lp;
	feasible = lhs_best >= rhs;
      default:
	fprintf(stderr, "Uncaught sense case %c\n", sense);
	rval = 1;
	goto CLEANUP;
      }

      epsilon_feasible = feasible && (fabs(lhs_best - rhs) < LP::EPSILON);
      printf("    CUT LHS -- best: %.2f, lp %.2f, rhs %.2f, sense %c, ",
	     lhs_best, lhs_lp, rhs, sense);
      printf("viol %.2f, nz %d\n", viol, num_nonzero);
      printf("    CUT IS... violated by LP solution:   %d\n", lp_violated);
      printf("              tight at best tour:        %d\n", exact_cut);
      printf("              feasible at best tour:     %d\n", feasible);
      printf("              eps-feasible at best tour: %d\n", epsilon_feasible);

      if(!feasible || !lp_violated){
	printf("    discarding infeasible or unviolated cut\n");
	continue;
      }
      
      if(!epsilon_feasible){
	printf("   discarding non-tight, non epsilon-feasible cut\n");
	continue;
      }

      if(arg->is_best_exact && !exact_cut){
	printf("   discarding bc exact cuts take precedence over inexact\n");
	continue;
      }

      
      if( (!(arg->is_best_exact) && exact_cut) ||
	  (viol > arg->best_viol)){
	arg->best_nonzeros = num_nonzero;
	arg->best_sense = sense;
	arg->best_viol = viol;
	arg->best_rhs = rhs;
	arg->is_best_exact = exact_cut;
	for(int k = 0; k < num_nonzero; k++){
	  arg->best_coeffs[k] = arg->coefficient_buffer[k];
	  arg->best_indices[k] = arg->index_buffer[k];
	}

	printf("    ****NEW best cut exact: %d, viol: %.2f\n", exact_cut, viol);
      }
    }
  }

  // if(arg->is_best_exact){
  //   printf("    found exact best cut...trying to prematurely terminate?\n");
  //   *useraction_p = CPX_CALLBACK_FAIL;
  // }
  
 CLEANUP:
  if(rval)
    fprintf(stderr, " problem in Cut<general>::solvecallback\n");
  return rval;
}
