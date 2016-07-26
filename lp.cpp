#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<getopt.h>
#include<cplex.h>
#include "lp.h"

int PSEPlp_init (PSEPlp *lp){
  int rval = 0;
  
  lp->cplex_lp = (CPXLPptr) NULL;
  lp->cplex_env = CPXopenCPLEX(&rval);

  if(rval){
    fprintf (stderr, "CPXopenCPLEX failed, return code %d\n", rval);
    goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

void PSEPlp_free (PSEPlp *lp){
  if (lp) {
    if (lp->cplex_env){
      if (lp->cplex_lp){
	CPXfreeprob (lp->cplex_env, &(lp->cplex_lp));
      }
    CPXcloseCPLEX (&lp->cplex_env);
    lp->cplex_env = (CPXENVptr) NULL;
    }
  }
}

int PSEPlp_create (PSEPlp *lp, const char *name){
  int rval = 0;
  char nambuf[32];
  if (!lp->cplex_env){
    fprintf (stderr, "no cplex_env in lp\n");
    rval = 1; goto CLEANUP;
  }

  lp->cplex_lp = CPXcreateprob (lp->cplex_env, &rval, nambuf);
  if (!lp->cplex_lp || rval){
    fprintf (stderr, "CPXcreateprob failed, return code %d\n", rval);
  }

 CLEANUP:
  return rval;
}

int PSEPlp_new_row (PSEPlp *lp, char sense, double rhs){
  int rval = 0;
  char asense[1];
  double arhs[1];

  asense[0] = sense;
  arhs[0] = rhs;

  rval = CPXnewrows (lp->cplex_env, lp->cplex_lp, 1, arhs, asense,
		     (double *) NULL, (char **) NULL);
  if (rval) {
    fprintf (stderr, "CPXnewrows failed\n"); goto CLEANUP;
  }
 CLEANUP:
  return rval;
}

int PSEPlp_addrows (PSEPlp *lp, int newrows, int newnz, double *rhs,
		    char *sense, int *rmatbeg, int *rmatind,
		    double *rmatval){
  int rval = 0;

  rval = CPXaddrows (lp->cplex_env, lp->cplex_lp, 0, newrows, newnz,
		     rhs, sense, rmatbeg, rmatind, rmatval,
		     (char **) NULL, (char **) NULL);
  if (rval) {fprintf(stderr, "CPXaddrows failed\n"); goto CLEANUP;}

 CLEANUP:
  return rval;
}

int PSEPlp_getrows(PSEPlp *lp, int *nzcnt_p, int *rmatbeg, int *rmatind,
		    double *rmatval, int rmatspace, int *surplus_p, int begin,
		   int end){
  int rval = CPXgetrows(lp->cplex_env, lp->cplex_lp, nzcnt_p, rmatbeg,
			rmatind, rmatval, rmatspace, surplus_p, begin, end);
  /*if(rval){
    if (rval == CPXERR_NEGATIVE_SURPLUS)
      fprintf(stderr, "Insufficient space in rmatind, rmatval ");
    fprintf(stderr, "PSEPlp_getrows failed\n");
    }*/

  return rval;
}

int PSEPlp_delrows(PSEPlp *lp, int begin, int end){
  int rval = CPXdelrows(lp->cplex_env, lp->cplex_lp, begin, end);
  if(rval) {fprintf(stderr, "CPXdelrows failed\n");}
  return rval;
}

int PSEPlp_delsetrows(PSEPlp *lp, int *delstat){
  int rval = CPXdelsetrows(lp->cplex_env, lp->cplex_lp, delstat);
  if(rval) {fprintf(stderr, "CPXdelsetrows failed, rval %d\n", rval);}
  return rval;
}

int PSEPlp_numrows(PSEPlp *lp){
  return CPXgetnumrows (lp->cplex_env, lp->cplex_lp);
}

int PSEPlp_addcols (PSEPlp *lp, int newcols, int newnz, double *obj,
		    int *cmatbeg, int *cmatind, double *cmatval,
		    double *lb, double *ub){
  int rval = 0;

  rval = CPXaddcols (lp->cplex_env, lp->cplex_lp, newcols, newnz, obj,
		     cmatbeg, cmatind, cmatval, lb, ub,
		     (char **) NULL);
  if (rval) {fprintf (stderr, "CPXaddcols failed\n"); goto CLEANUP;}

 CLEANUP:
  return rval;
}

int PSEPlp_setbnd (PSEPlp *lp, int col, char lower_or_upper,
		   double bnd){
  int rval = 0;
  int cindex[1];
  double bd[1];
  char lu[1];
  cindex[0] = col;
  lu[0] = lower_or_upper;
  bd[0] = bnd;

  rval = CPXchgbds (lp->cplex_env, lp->cplex_lp, 1, cindex, lu, bd);
  
  if (rval){
    fprintf (stderr, "couldn't set bnd on var %d in cplex\n", col);
    goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_clampbnd(PSEPlp *lp, int col, char lower_or_upper, double bnd){
  int rval = 0, count = 1;
  
  rval = CPXtightenbds(lp->cplex_env, lp->cplex_lp, count, &col,
		       &lower_or_upper, &bnd);
  if(rval)
    fprintf(stderr, "Problem clamping with CPXtightenbds, rval %d\n", rval);
  return rval;
}

int PSEPlp_relaxbds(PSEPlp *lp, int count, int const *indices,
		    char const *lower_or_upper, double const *bd){
  int rval = CPXtightenbds(lp->cplex_env, lp->cplex_lp, count, indices,
			   lower_or_upper, bd);
  if(rval)
    fprintf(stderr, "Problem relaxing w CPXtightenbds, rval %d\n", rval);
  return rval;
}

int PSEPlp_no_opt (PSEPlp *lp){
  int rval = 0, solstat;

  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_Simplex_Limits_Iterations, 0);
  if (rval){
    fprintf (stderr, "Failed to set limit to 0\n");
    goto CLEANUP;
  }

  rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
  if (rval){
    fprintf(stderr, "CPXdualopt zero iteration failed\n");
    goto CLEANUP;
  }

  solstat = CPXgetstat(lp->cplex_env, lp->cplex_lp);
  if (solstat != CPX_STAT_ABORT_IT_LIM){
    fprintf(stderr, "Failed to perform zero iterations\n");
    rval = 1;
    goto CLEANUP;
  } /*else {
    printf("Successfully performed zero iterations\n");
    }*/
  
  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_Simplex_Limits_Iterations,
			 LP::DEFAULT_ITLIM);
  if (rval){
    fprintf (stderr, "Failed to revert limit\n");
    goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_primal_opt (PSEPlp *lp, int *infeasible){
  int rval = 0, solstat;

  *infeasible = 0;

  rval = CPXprimopt (lp->cplex_env, lp->cplex_lp);
  if (rval){
    fprintf (stderr, "CPXprimopt failed\n");
    goto CLEANUP;
  }

  solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
  if (solstat == CPX_STAT_INFEASIBLE) {
    if (infeasible) *infeasible = 1;
  }/* else if (solstat == CPX_STAT_ABORT_IT_LIM){
    printf("Successfully performed one iteration of simplex\n");
    }*/ else if (solstat != CPX_STAT_OPTIMAL &&
	     solstat != CPX_STAT_OPTIMAL_INFEAS) {
    fprintf (stderr, "Cplex opt status %d\n", solstat);
    rval = 1; goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_dual_opt (PSEPlp *lp, int *infeasible){
  int rval = 0, solstat;

  *infeasible = 0;

  rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
  if (rval){
    fprintf (stderr, "CPXdualopt failed\n");
    goto CLEANUP;
  }

  solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
  if (solstat == CPX_STAT_INFEASIBLE) {
    if(infeasible)
      *infeasible = 1;
    fprintf(stderr, "PSEPlp_dual_opt detected infeasibility\n");
    goto CLEANUP;
  }else if (solstat != CPX_STAT_OPTIMAL &&
	     solstat != CPX_STAT_OPTIMAL_INFEAS) {
    fprintf (stderr, "Cplex opt status %d\n", solstat);
    rval = 1; goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_primal_pivot (PSEPlp *lp, int *infeasible){
  int rval = 0, solstat;

  
  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_Simplex_Limits_Iterations, 1);
  if (rval){
    fprintf (stderr, "Failed to set limit to 1\n");
    goto CLEANUP;
  }
  
  rval = CPXprimopt (lp->cplex_env, lp->cplex_lp);
  if (rval){
    fprintf (stderr, "CPXprimopt failed\n");
    goto CLEANUP;
  }

  solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
  if (solstat == CPX_STAT_INFEASIBLE){
    if (infeasible)
      *infeasible = 1;
  } else if (solstat != CPX_STAT_OPTIMAL &&
	     solstat != CPX_STAT_OPTIMAL_INFEAS &&
	     solstat != CPX_STAT_ABORT_IT_LIM){
    fprintf (stderr, "Cplex opt satus %d\n", solstat);
    rval = 1;
    goto CLEANUP;
  }
  
  if (rval){
    fprintf (stderr, "Failed to perform simplex pivot\n");
    goto CLEANUP;
  }

  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_Simplex_Limits_Iterations,
			 LP::DEFAULT_ITLIM);
  if (rval){
    fprintf (stderr, "Failed to revert itlimit\n");
    goto CLEANUP;
  }
  


 CLEANUP:
  return rval;
}

int PSEPlp_dual_pivot (PSEPlp *lp, int *infeasible){
  int rval = 0, solstat;

  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_Simplex_Limits_Iterations, 1);
  if (rval){
    fprintf (stderr, "Failed to set limit to 1\n");
    goto CLEANUP;
  }
  
  rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
  if (rval){
    fprintf (stderr, "CPXprimopt failed\n");
    goto CLEANUP;
  }

  solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
  if (solstat == CPX_STAT_INFEASIBLE){
    if (infeasible)
      *infeasible = 1;
  } else if (solstat != CPX_STAT_OPTIMAL &&
	     solstat != CPX_STAT_OPTIMAL_INFEAS &&
	     solstat != CPX_STAT_ABORT_IT_LIM){
    fprintf (stderr, "Cplex opt satus %d\n", solstat);
    rval = 1;
    goto CLEANUP;
  }
  
  if (rval){
    fprintf (stderr, "Failed to perform simplex pivot\n");
    goto CLEANUP;
  }

  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_Simplex_Limits_Iterations,
			 LP::DEFAULT_ITLIM);
  if (rval){
    fprintf (stderr, "Failed to revert itlimit\n");
    goto CLEANUP;
  }


 CLEANUP:
  return rval;
}

int PSEPlp_pivot (PSEPlp *lp, const int entering_var, const int leaving_var,
		  const int nb_status){
  int rval = CPXpivot (lp->cplex_env, lp->cplex_lp, entering_var, leaving_var,
		       nb_status);
  if(rval)
    fprintf (stderr, "CPXpivot failed, rval %d\n", rval);

  return rval;
}

int PSEPlp_getobj (PSEPlp *lp, double *obj, int numcols){
  int rval = CPXgetobj(lp->cplex_env, lp->cplex_lp, obj, 0, numcols - 1);
  if(rval)
    fprintf(stderr, "CPXgetobj failed, rval %d\n", rval);
  return rval;
}

int PSEPlp_objval (PSEPlp *lp, double *obj){
  int rval = 0;

  rval = CPXgetobjval (lp->cplex_env, lp->cplex_lp, obj);
  if (rval) {
    fprintf (stderr, "CPXgetobjval failed \n"); goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_x (PSEPlp *lp, double *x){
  int rval = 0, ncols;
  ncols = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
  if (ncols == 0){
    fprintf (stderr, "No cols in LP\n");
    rval = 1; goto CLEANUP;
  }
  rval = CPXgetx (lp->cplex_env, lp->cplex_lp, x, 0, ncols - 1);
  if (rval) {
    fprintf (stderr, "CPXgetx failed\n"); goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_write (PSEPlp *lp, const char *fname){
  int rval = 0;
  char nambuf[32], lpbuf[4];

  strncpy (nambuf, fname, sizeof(nambuf));
  nambuf[sizeof(nambuf) - 1] = '\0';
  strncpy (lpbuf, "RLP", sizeof ("RLP"));

  rval = CPXwriteprob (lp->cplex_env, lp->cplex_lp, nambuf, lpbuf);
  if (rval){
    fprintf (stderr, "CPXwriteprob failed\n"); goto CLEANUP;
  }

 CLEANUP:
  return rval;
}

int PSEPlp_dualfeas(PSEPlp *lp){
  int dualfeas, rval = 0, result = 0;

  rval = CPXsolninfo (lp->cplex_env, lp->cplex_lp, NULL, NULL, NULL,
		      &dualfeas);

  if (rval){
    fprintf (stderr, "CPXsolninfo failed\n");
    return 0;
  }

  result = dualfeas;
  return result;
}

int PSEPlp_solstat(PSEPlp *lp){
  return CPXgetstat(lp->cplex_env, lp->cplex_lp);
}

int PSEPlp_chgobj (PSEPlp *lp, int count, int const * indices,
		  double const * values){
  int rval = 0;
  rval = CPXchgobj(lp->cplex_env, lp->cplex_lp, count, indices, values);
  if(rval){
    fprintf (stderr, "chgobj failed\n");
    return rval;
  } else
    return rval;
}

int PSEPlp_copybase (PSEPlp *lp, int *colstat, int *rowstat){
  int rval = 0;
  rval = CPXcopybase(lp->cplex_env, lp->cplex_lp, colstat, rowstat);
  if(rval){
    fprintf (stderr, "copybase failed, rval %d\n", rval);
  }
  return rval;
}

int PSEPlp_getbase (PSEPlp *lp, int * colstat, int * rowstat){
  int rval = 0;
  if(!colstat && !rowstat){
    fprintf(stderr, "Passed two null arrays to getbase\n");
    return rval;
  }
  
  rval = CPXgetbase(lp->cplex_env, lp->cplex_lp, colstat, rowstat);
  if (rval){
    fprintf(stderr, "getbase failed, rval %d\n", rval);
  }

  return rval;
}

int PSEPlp_getslack(PSEPlp *lp, double *slack, int begin, int end){
  int rval = 0;
  rval = CPXgetslack(lp->cplex_env, lp->cplex_lp, slack, begin, end);
  if(rval)
    fprintf(stderr, "lp_getslack failed, rval %d\n", rval);

  return rval;
}

int PSEPlp_bhead (PSEPlp *lp, int *head, double *x){
  int rval = CPXgetbhead(lp->cplex_env, lp->cplex_lp, head, x);

  if(rval)
    fprintf(stderr, "bhead failed, rval %d\n", rval);

  return rval;
}

int PSEPlp_getsense (PSEPlp *lp, char *sense, int rownum){
  int rval = CPXgetsense(lp->cplex_env, lp->cplex_lp, sense, rownum, rownum);

  if(rval)
    fprintf(stderr, "getsense failed, rval %d\n", rval);

  return rval;
}

int PSEPlp_get_redcosts (PSEPlp *lp, double * cost_array){
  int rval = 0;
  int numcols = CPXgetnumcols(lp->cplex_env, lp->cplex_lp);
  rval = CPXgetdj(lp->cplex_env, lp->cplex_lp, cost_array, 0, numcols - 1);

  if(rval){
    fprintf(stderr, "get_redcosts failed\n");
  }

  return rval;
}
