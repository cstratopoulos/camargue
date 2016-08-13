#include<stdio.h>

#include "mip.h"

int PSEPlp_make_mip (PSEPlp *lp){
  int rval = CPXchgprobtype(lp->cplex_env, lp->cplex_lp, CPXPROB_MILP);
  if(rval)
    fprintf(stderr, "PSEPlp_make_mip failed, rval %d\n", rval);
  return rval;
}

int PSEPlp_make_lp (PSEPlp *lp){
    int rval = CPXchgprobtype(lp->cplex_env, lp->cplex_lp, CPXPROB_LP);
  if(rval)
    fprintf(stderr, "PSEPlp_make_lp failed, rval %d\n", rval);
  return rval;
}

int PSEPmip_param (PSEPlp *lp){
  int rval = 0;

  //MIP PARAMS
  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Limits_Nodes, 0);
  if(rval){ fprintf(stderr, "MIP nodelimit"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Strategy_Search,
			CPX_MIPSEARCH_TRADITIONAL);
  if(rval){ fprintf(stderr, "MIPSEARCH"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Strategy_StartAlgorithm,
			CPX_ALG_PRIMAL);
  if(rval){ fprintf(stderr, "MIP root algorithm"); goto CLEANUP; }

  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_MIP_Strategy_RINSHeur, -1);
  if(rval){ fprintf(stderr, "RINSEHEUR"); goto CLEANUP; }
  
  rval = CPXsetlongparam(lp->cplex_env, CPXPARAM_MIP_Strategy_HeuristicFreq,
			 -1);
  if(rval){ fprintf(stderr, "HEURFREQ"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Strategy_FPHeur, -1);
  if(rval){ fprintf(stderr, "FPHEUR"); goto CLEANUP; }
  
  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Strategy_Probe, -1);
  if(rval){ fprintf(stderr, "PROBE"); goto CLEANUP; }


  //PREPROCESSING/PRESOLVE
  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_Preprocessing_Presolve, 0);
  if(rval){ fprintf(stderr, "Presolve"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_Preprocessing_Relax, 0);
  if(rval){ fprintf(stderr, "Presolve LP relaxation"); goto CLEANUP; }

  //CUTS
  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_Cliques, -1);
  if(rval){ fprintf(stderr, "Clique cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_Covers, -1);
  if(rval){ fprintf(stderr, "Cover cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_Disjunctive, -1);
  if(rval){ fprintf(stderr, "Disjunctive cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_FlowCovers, -1);
  if(rval){ fprintf(stderr, "Flow cover cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_PathCut, -1);
  if(rval){ fprintf(stderr, "Flow path cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_Gomory, -1);
  if(rval){ fprintf(stderr, "Gomory frac cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_GUBCovers, -1);
  if(rval){ fprintf(stderr, "GUB cover cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_Implied, -1);
  if(rval){ fprintf(stderr, "Implied cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_LocalImplied, -1);
  if(rval){ fprintf(stderr, "Local Implied cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_LiftProj, -1);
  if(rval){ fprintf(stderr, "Lift and Project cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_MCFCut, -1);
  if(rval){ fprintf(stderr, "MCF cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_MIRCut, -1);
  if(rval){ fprintf(stderr, "MIR cuts"); goto CLEANUP; }

  rval = CPXsetintparam(lp->cplex_env, CPXPARAM_MIP_Cuts_ZeroHalfCut, -1);
  if(rval){ fprintf(stderr, "0-1/2 cuts"); goto CLEANUP; }
  

 CLEANUP:
  if(rval)
    fprintf(stderr, " failed in PSEPlp_mip_param, rval %d\n", rval);

  return rval;
}

int PSEPmip_change_vartype (PSEPlp *lp, int count, int const *indices,
			char const *vartype){
  int rval = CPXchgctype(lp->cplex_env, lp->cplex_lp, count, indices, vartype);
  if(rval)
    fprintf(stderr, "PSEPlp_change_vartype failed, rval %d\n", rval);
  return rval;
}

int PSEPmip_opt (PSEPlp *lp){
  int rval = CPXmipopt(lp->cplex_env, lp->cplex_lp);
  if(rval) fprintf(stderr, "PSEPlp_mip_opt failed, rval %d\n", rval);
  return rval;
}

int PSEPmip_getbestobjval (PSEPlp *lp, double *obj){
  int rval = CPXgetbestobjval(lp->cplex_env, lp->cplex_lp, obj);
  if(rval)
    fprintf(stderr, "PSEPlp_getbestobjval failed\n");
  return rval;
}

int PSEPmip_getnumcuts (PSEPlp *lp, int cuttype, int *num_p){
  int rval = CPXgetnumcuts(lp->cplex_env, lp->cplex_lp, cuttype, num_p);
  if(rval)
    fprintf(stderr, "PSEPlp_getnumcuts failed, rval %d\n", rval);
  return rval;
}
