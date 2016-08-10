#ifndef PSEP_MIP_H
#define PSEP_MIP_H

#include "lp.h"

//change LP to MILP problem type
int PSEPlp_make_mip (PSEPlp *lp);

//change MIP back to LP problem type
int PSEPlp_make_lp (PSEPlp *lp);

//disables all of the CPLEX heuristics, all cuts except disjunctive/frac/GMI,
//sets node limit to zero
int PSEPmip_param (PSEPlp *lp);

//makes all the variables binary in a MIP
int PSEPmip_change_vartype (PSEPlp *lp, int count, int const *indices,
			char const *vartype);

//calls CPXmipopt on a MIP problem
int PSEPmip_opt (PSEPlp *lp);

//get best objval of a mip object
int PSEPmip_getbestobjval (PSEPlp *lp, double *obj);

//calls CPXgetnumcuts to get the number of gomory fractional cuts
int PSEPmip_getnumcuts (PSEPlp *lp, int cuttype, int *num_p);

#endif
