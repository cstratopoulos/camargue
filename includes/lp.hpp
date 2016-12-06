#ifndef CMR_LP_H
#define CMR_LP_H

#include <cplex.h>
#include "util.hpp"
                
//LP object consisting of CPLEX environment pointer and CPLEX LP pointer
struct CMRlp {
  CPXENVptr cplex_env;
  CPXLPptr cplex_lp;
};

//all functions return an int rval that acts as a status message--an rval of
//0 indicates successful function call, other values signify errors. 


//interface function to open cplex and create lp environment
int CMRlp_init (CMRlp *lp);

//interface to CPXcreateprob--creates an empty LP problem in the environment
int CMRlp_create (CMRlp *lp, const char *name);

//write LP to file
int CMRlp_write (CMRlp *lp, const char *fname);

//delete/free an LP environment
void CMRlp_free (CMRlp *lp);

//creates a new row with empty constraints in lp problem object
// sense L E G for <=, =, >=, rhs is righthand-side of constraint
int CMRlp_new_row (CMRlp *lp, char sense, double rhs);

//adds rows to lp problem object, sense, rhs, as above
//newrows is # of new rows
//rmatbeg/rmatind are used with rmatval to define rows to be added
//newnz is number of new nonzero constraints
int CMRlp_addrows (CMRlp *lp, int newrows, int newnz, double *rhs,
		    char *sense, int *rmatbeg, int *rmatind,
		    double *rmatval);

int CMRlp_getrows (CMRlp *lp, int *nzcnt_p, int *rmatbeg, int *rmatind,
		    double *rmatval, int rmatspace, int *surplus_p, int begin,
		    int end);

//returns the number of rows/cols in lp problem object
int CMRlp_numrows(CMRlp *lp);
int CMRlp_numcols(CMRlp *lp);


//deletes specified range of rows
int CMRlp_delrows(CMRlp *lp, int begin, int end);

//deletes set of rows (resp cols) indicated by delstat, deleting those where
//delstat[i] = 1
//afterwards delstat[i] is
// -- -1 if the row (resp col) was deleted
// -- the new index number assigned if not
int CMRlp_delsetrows(CMRlp *lp, int *delstat);
int CMRlp_delsetcols(CMRlp *lp, int *delstat);

//similarly as above
//obj is an array of length newcols containing objective function coeffs
//for new variables
//lb and ub are arrays of length newcols containing lower and upper bounds
//on the new vars
int CMRlp_addcols (CMRlp *lp, int newcols, int newnz, double *obj,
		    int *cmatbeg, int *cmatind, double *cmatval,
		    double *lb, double *ub);

//used to change lower_or_upper bound on one column with index col of the LP
//lower_or_upper is L, U, or B for upper and lower. 
int CMRlp_setbnd (CMRlp *lp, int col, char lower_or_upper,
		   double bnd);

int CMRlp_clampbnd (CMRlp *lp, int col, char lower_or_upper, double bnd);
int CMRlp_relaxbds (CMRlp *lp, int count, int const *indices,
		     char const *lower_or_upper, double const * bd);


//performs zero primal simplex pivots but factors the basis
int CMRlp_no_opt (CMRlp *lp);

//finds a solution to lp by calling cplex dual simplex algorithm
//infeasible used to return whether the lp was found infeasible
int CMRlp_dual_opt (CMRlp *lp, int *infeasible);
int CMRlp_primal_opt (CMRlp *lp, int *infeasible);

//perform a single primal simplex pivot
int CMRlp_primal_pivot (CMRlp *lp, int *infeasible);

//performs non-degenerate primal simplex pivot
int CMRlp_primal_nd_pivot (CMRlp *lp, int *infeasible, const double lowlimit);

//perform a single dual pivot of the simplex algorithm
int CMRlp_dual_pivot (CMRlp *lp, int *infeasible);

//gets the objective function
int CMRlp_getobj (CMRlp *lp, double *obj, int numcols);

//gets number of simplex iterations
int CMRlp_itcount (CMRlp *lp);

//calls CPXgetobjval to obtain the optimal objective value of lp, with a
//pointer to obj, where the result is to be stored
int CMRlp_objval (CMRlp *lp, double *obj);

//argmin of objval, if successful stores the optimal solution in an array
//x of type double
int CMRlp_x (CMRlp *lp, double *x);

//bool to return if solution is dual feasible
int CMRlp_dualfeas (CMRlp *lp);

int CMRlp_solstat (CMRlp *lp);

//changes the objective function of existing LP, used to feed dummy solution
//and then revert to actual LP
int CMRlp_chgobj (CMRlp *lp, int count, int const * indices,
		  double const * values);

//change the sense of the constraints to sense
int CMRlp_chgsense (CMRlp *lp, const int count, int const * indices,
		     char const * sense);

//change an individual coefficient, col = -1 for RHS
int CMRlp_chgcoef (CMRlp *lp, const int row, const int col,
		    const double newvalue);

//copies a starting solution and/or basis for use by optimizers
int CMRlp_copystart (CMRlp *lp, int const * cstat, int const * rstat,
		      double const * cprim, double const * rprim,
		      double const * cdual, double const * rdual);

//copies basis statuses into the LP
int CMRlp_copybase ( CMRlp *lp, int *colstat, int *rowstat);

//access resident basis. one of colstat or rowstat may be NULL if not needed
int CMRlp_getbase (CMRlp *lp, int * colstat, int *rowstat);

//access the slacks for a range of rows
int CMRlp_getslack (CMRlp *lp, double *slack, int begin, int end);

//access the dual variables for a range of constraints
int CMRlp_getpi (CMRlp *lp, double *pi, int begin, int end);

//access the lower bounds on a range of variables
int CMRlp_getlb (CMRlp *lp, double *lb, int begin, int end);

//computes infeasibility of a given solution
//the array feas_stat will be nonzero in entry i if the ith row is violated
//by x. If x is null, the current LP solution is used.
int CMRlp_getrowinfeas (CMRlp *lp, double const *x, double *feas_stat,
			 int begin, int end);

//access the basis header, stored in head, and basic variable values stored in
//x; either may be NULL if not needed
int CMRlp_bhead (CMRlp *lp, int *head, double *x);

//access the sense of the constraint in rownum:
// 'L' for <=
// 'G' for >=
int CMRlp_getsense (CMRlp *lp, char *sense, int rownum);

//access array of reduced costs
int CMRlp_get_redcosts (CMRlp *lp, double * cost_array);

#endif
