#ifndef PSEP_LP_H
#define PSEP_LP_H

#include <cplex.h>
#include "PSEP_util.h"
                
//LP object consisting of CPLEX environment pointer and CPLEX LP pointer
struct PSEPlp {
  CPXENVptr cplex_env;
  CPXLPptr cplex_lp;
};

//all functions return an int rval that acts as a status message--an rval of
//0 indicates successful function call, other values signify errors. 


//interface function to open cplex and create lp environment
int PSEPlp_init (PSEPlp *lp);

//interface to CPXcreateprob--creates an empty LP problem in the environment
int PSEPlp_create (PSEPlp *lp, const char *name);

//write LP to file
int PSEPlp_write (PSEPlp *lp, const char *fname);

//delete/free an LP environment
void PSEPlp_free (PSEPlp *lp);

//creates a new row with empty constraints in lp problem object
// sense L E G for <=, =, >=, rhs is righthand-side of constraint
int PSEPlp_new_row (PSEPlp *lp, char sense, double rhs);

//adds rows to lp problem object, sense, rhs, as above
//newrows is # of new rows
//rmatbeg/rmatind are used with rmatval to define rows to be added
//newnz is number of new nonzero constraints
int PSEPlp_addrows (PSEPlp *lp, int newrows, int newnz, double *rhs,
		    char *sense, int *rmatbeg, int *rmatind,
		    double *rmatval);

int PSEPlp_getrows (PSEPlp *lp, int *nzcnt_p, int *rmatbeg, int *rmatind,
		    double *rmatval, int rmatspace, int *surplus_p, int begin,
		    int end);

//returns the number of rows/cols in lp problem object
int PSEPlp_numrows(PSEPlp *lp);
int PSEPlp_numcols(PSEPlp *lp);


//deletes specified range of rows
int PSEPlp_delrows(PSEPlp *lp, int begin, int end);

//deletes set of rows (resp cols) indicated by delstat, deleting those where
//delstat[i] = 1
//afterwards delstat[i] is
// -- -1 if the row (resp col) was deleted
// -- the new index number assigned if not
int PSEPlp_delsetrows(PSEPlp *lp, int *delstat);
int PSEPlp_delsetcols(PSEPlp *lp, int *delstat);

//similarly as above
//obj is an array of length newcols containing objective function coeffs
//for new variables
//lb and ub are arrays of length newcols containing lower and upper bounds
//on the new vars
int PSEPlp_addcols (PSEPlp *lp, int newcols, int newnz, double *obj,
		    int *cmatbeg, int *cmatind, double *cmatval,
		    double *lb, double *ub);

//used to change lower_or_upper bound on one column with index col of the LP
//lower_or_upper is L, U, or B for upper and lower. 
int PSEPlp_setbnd (PSEPlp *lp, int col, char lower_or_upper,
		   double bnd);

int PSEPlp_clampbnd (PSEPlp *lp, int col, char lower_or_upper, double bnd);
int PSEPlp_relaxbds (PSEPlp *lp, int count, int const *indices,
		     char const *lower_or_upper, double const * bd);


//performs zero primal simplex pivots but factors the basis
int PSEPlp_no_opt (PSEPlp *lp);

//finds a solution to lp by calling cplex dual simplex algorithm
//infeasible used to return whether the lp was found infeasible
int PSEPlp_dual_opt (PSEPlp *lp, int *infeasible);
int PSEPlp_primal_opt (PSEPlp *lp, int *infeasible);

//perform a single primal simplex pivot
int PSEPlp_primal_pivot (PSEPlp *lp, int *infeasible);

//performs non-degenerate primal simplex pivot
int PSEPlp_primal_nd_pivot (PSEPlp *lp, int *infeasible, const double lowlimit);

//perform a single dual pivot of the simplex algorithm
int PSEPlp_dual_pivot (PSEPlp *lp, int *infeasible);

//gets the objective function
int PSEPlp_getobj (PSEPlp *lp, double *obj, int numcols);

//gets number of simplex iterations
int PSEPlp_itcount (PSEPlp *lp);

//calls CPXgetobjval to obtain the optimal objective value of lp, with a
//pointer to obj, where the result is to be stored
int PSEPlp_objval (PSEPlp *lp, double *obj);

//argmin of objval, if successful stores the optimal solution in an array
//x of type double
int PSEPlp_x (PSEPlp *lp, double *x);

//bool to return if solution is dual feasible
int PSEPlp_dualfeas (PSEPlp *lp);

int PSEPlp_solstat (PSEPlp *lp);

//changes the objective function of existing LP, used to feed dummy solution
//and then revert to actual LP
int PSEPlp_chgobj (PSEPlp *lp, int count, int const * indices,
		  double const * values);

//change the sense of the constraints to sense
int PSEPlp_chgsense (PSEPlp *lp, const int count, int const * indices,
		     char const * sense);

//change an individual coefficient, col = -1 for RHS
int PSEPlp_chgcoef (PSEPlp *lp, const int row, const int col,
		    const double newvalue);

//copies a starting solution and/or basis for use by optimizers
int PSEPlp_copystart (PSEPlp *lp, int const * cstat, int const * rstat,
		      double const * cprim, double const * rprim,
		      double const * cdual, double const * rdual);

//copies basis statuses into the LP
int PSEPlp_copybase ( PSEPlp *lp, int *colstat, int *rowstat);

//access resident basis. one of colstat or rowstat may be NULL if not needed
int PSEPlp_getbase (PSEPlp *lp, int * colstat, int *rowstat);

//access the slacks for a range of rows
int PSEPlp_getslack (PSEPlp *lp, double *slack, int begin, int end);

//access the dual variables for a range of constraints
int PSEPlp_getpi (PSEPlp *lp, double *pi, int begin, int end);

//access the lower bounds on a range of variables
int PSEPlp_getlb (PSEPlp *lp, double *lb, int begin, int end);

//access the basis header, stored in head, and basic variable values stored in
//x; either may be NULL if not needed
int PSEPlp_bhead (PSEPlp *lp, int *head, double *x);

//access the sense of the constraint in rownum:
// 'L' for <=
// 'G' for >=
int PSEPlp_getsense (PSEPlp *lp, char *sense, int rownum);

//access array of reduced costs
int PSEPlp_get_redcosts (PSEPlp *lp, double * cost_array);

#endif
