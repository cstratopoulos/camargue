#include<algorithm>
#include<iostream>
#include<iomanip>
#include<tuple>
#include<utility>

#include<cmath>

#define DO_SAFE_MIR_DBL 1
#define SAFE_MIR_DEBUG_LEVEL DBG_LEVEL_HIGH
#define CUTSslackSign( row ) ( row->sense == 'L' ? 1 : (row->sense == 'E' ? 0 : -1 ) )
#include<safemir/src/gen_slvr.hpp>
#include<safemir/src/util_cuts.hpp>
#include<safemir/src/gen_mir.cpp>
#include<safemir/src/cplex_slvr.cpp>
#include<safemir/src/ds_slvr.cpp>
#include<safemir/src/ds_cuts.cpp>
#include<safemir/src/safe_mir_dbl.cpp>

#include "safegmi.hpp"

using namespace std;
using namespace PSEP;

int Cut<safeGMI>::cutcall(){
  int rval = 0;
  
  rval = init_constraint_info();
  if(rval) goto CLEANUP;

  rval = get_tab_rows();
  if(rval) goto CLEANUP;

  rval = separate();
  if(rval) goto CLEANUP;

  rval = add_cut();
  if(rval) goto CLEANUP;
  
 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<safeGMI>::cutcall\n";
  safe_mir_data.reset(NULL);
  return rval;
}

int Cut<safeGMI>::init_constraint_info(){
  int rval = 0;

  // This should no longer be necessary with new cut add control flow
  // rval = PSEPlp_copystart(&m_lp, &frac_colstat[0], &frac_rowstat[0],
  // 			  &m_lp_edges[0], NULL, NULL, NULL);
  // if(rval) PSEP_GOTO_CLEANUP("Failed to copy frac solution, ");

  // rval = PSEPlp_no_opt(&m_lp);
  // if(rval) PSEP_GOTO_CLEANUP("Failed to factor basis, ");
  
  try { safe_mir_data.reset(new SafeMIRGroup(m_lp)); }
  catch(const std::bad_alloc &e){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for safe_mir_data reset, ");
  }
  
  rval = SLVRformulationRows(&(safe_mir_data->lp_obj),
			     &(safe_mir_data->constraint_matrix));
  if(rval) PSEP_GOTO_CLEANUP("SLVRformulationRows failed, ");

  rval = SLVRgetBasisInfo(&(safe_mir_data->lp_obj),
			  &(safe_mir_data->basis_info));
  if(rval) PSEP_GOTO_CLEANUP("SLVRgetBasisInfo failed, ");

  rval = SLVRgetVarInfo(&(safe_mir_data->lp_obj), true,
			&(safe_mir_data->var_info));
  if(rval) PSEP_GOTO_CLEANUP("SLVRgetVarInfo failed, ");

  safe_mir_data->full_x = SLVRgetFullX(&(safe_mir_data->lp_obj),
				       safe_mir_data->constraint_matrix,
				       &m_lp_edges[0]);
  if(!safe_mir_data->full_x) PSEP_GOTO_CLEANUP("SLVRgetFullX failed, ");
  

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<safeGMI>::get_constraint_matrix\n";
  return rval;
}

int Cut<safeGMI>::get_tab_rows(){
  int rval = 0;
  int numcols = m_lp_edges.size();
  int numrows = PSEPlp_numrows(&m_lp);

  vector<double> &lp_vranking = safe_mir_data->lp_vranking;
  vector<pair<int, double>> frac_basic_vars;

  try{ lp_vranking.resize(numrows + numcols, -1.0); }
  catch(const std::bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for lp_vranking, ");
  }

  for(int i = 0; i < support_indices.size(); i++){
    int ind = support_indices[i];
    double lp_entry = m_lp_edges[ind];
    if(frac_colstat[ind] == CPX_BASIC && lp_entry < 1 - LP::EPSILON){
      lp_vranking[ind] = (-(lp_entry - 0.5) * (lp_entry - 0.5)) + 0.25;
      try{ frac_basic_vars.push_back(pair<int, double>(ind, lp_entry)); }
      catch(const std::bad_alloc &){
	rval = 1; PSEP_GOTO_CLEANUP("Out of memory for frac_basic_vars, ");
      }
    }
  }
  //fractional basic vars ranked by how close they are to 0.5
  sort(frac_basic_vars.begin(), frac_basic_vars.end(),
       [](const pair<int, double> &a, const pair<int, double> &b) -> bool {
	 return fabs(0.5 - a.second) < fabs(0.5 - b.second);
       });
  
  rval = CUTSnewSystem(&(safe_mir_data->tableau_rows),
		       frac_basic_vars.size());
  if(rval) PSEP_GOTO_CLEANUP("Out of memory for tableau rows, ");

  {
    CUTSsystem_t<double> *tab_system = safe_mir_data->tableau_rows;
    for(int i = 0; i < frac_basic_vars.size(); i++){
      int ind = frac_basic_vars[i].first;
      rval = SLVRgetTableauRow(&(safe_mir_data->lp_obj),
			       &(safe_mir_data->constraint_matrix),
			       &(tab_system->rows[tab_system->sys_rows]),
			       &(safe_mir_data->basis_info),
			       ind);
      if(rval) PSEP_GOTO_CLEANUP("SLVRgetTableauRow failed, ");

      tab_system->sys_rows += 1;
    }
  }
  

 CLEANUP:
  if(rval)
    cerr << "problem in Cut<safeGMI>::get_tab_rows\n";
  return rval;  
}

int Cut<safeGMI>::separate(){
  int rval = 0;

  int num_added = 0, total_num_added = 0;
  int numcols = PSEPlp_numcols(&m_lp);
  

  rval = CUTSnewRowList(&(safe_mir_data->generated_cuts));
  if(rval) PSEP_GOTO_CLEANUP("CUTSnewRowList failed, ");

  rval = SLVRcutter_iter(0, &(safe_mir_data->settings),
			 safe_mir_data->constraint_matrix,
			 safe_mir_data->tableau_rows,
			 &total_num_added, &num_added,
			 safe_mir_data->full_x,
			 safe_mir_data->var_info,
			 false, NULL,
			 numcols,
			 safe_mir_data->generated_cuts,
			 &safe_mir_data->lp_vranking[0]);
  if(rval) PSEP_GOTO_CLEANUP("SLVRcutter_iter failed, ");

  if(safe_mir_data->generated_cuts->size == 0){
    cout << "No safe mir cuts found\n";
    rval = 2;
  }

 CLEANUP:
  if(rval == 1)
    cerr << "problem in Cut<safeGMI>::separate\n";
  return rval;
}

int Cut<safeGMI>::add_cut(){
  int rval = 0;
  int num_added = 0;
  
  int numcols = PSEPlp_numcols(&m_lp), rmatbeg = 0;
  vector<double> best_edges;



  cut_obj best_cut(false, 0, -1.0, (CUTSsprow_t<double> *) NULL);
  
  try{
    for(int i = 0; i < best_tour_edges.size(); i++)
      best_edges.push_back(best_tour_edges[i]);
  } catch(const std::bad_alloc &){
    rval = 1; PSEP_GOTO_CLEANUP("Out of memory for integral tour vector, ")
  }

  for(CUTSrowListElem_t<double> *it = safe_mir_data->generated_cuts->first;
      it; it = it->next){
    CUTSsprow_t<double> *cur_row = it->row;

    int nz = cur_row->nz;
    double lp_viol, tour_act;
    double rhs = cur_row->rhs;
    if(cur_row->sense != 'G'){
      rval = 1; PSEP_GOTO_CLEANUP("Non geq cut??" );
    }
      
    rval = (CUTScomputeViolation(cur_row, &m_lp_edges[0], &lp_viol) ||
	    CUTScomputeActivity(cur_row, &best_edges[0], &tour_act));
    if(rval) PSEP_GOTO_CLEANUP("CUTScomputeActivity failed, ");

    bool exact = (tour_act == rhs);

    if((tour_act < rhs) || fabs(tour_act - rhs) >= LP::EPSILON)
      continue;

    cut_obj current_cut(exact, numcols - nz, lp_viol, cur_row);
    if(current_cut > best_cut){
      best_cut = current_cut;
      local_q.push_front(current_cut);
    }
  }
  

  if(!best_cut.row || local_q.empty()){
    rval = 2;
    cout << "No feasible cuts found\n";
    goto CLEANUP;
  }


  best_cut = local_q.peek_front();
  if(!best_cut.exact){
    rval = PSEPlp_addrows(&m_lp, 1, numcols - best_cut.zeros,
			  &best_cut.row->rhs, &best_cut.row->sense,
			  &rmatbeg, best_cut.row->rowind,
			  best_cut.row->rowval);
    PSEP_CHECK_RVAL(rval, "Couldn't add cut, ");

    num_added++;

    //TODO: replace this with a propper .clear() method
    while(!local_q.empty()) local_q.pop_front();
  }

  while(!local_q.empty()){
    best_cut = local_q.peek_front();

    if(best_cut.exact){
      rval = PSEPlp_addrows(&m_lp, 1, numcols - best_cut.zeros,
			    &best_cut.row->rhs, &best_cut.row->sense,
			    &rmatbeg, best_cut.row->rowind,
			    best_cut.row->rowval);
      PSEP_CHECK_RVAL(rval, "Couldn't add cut, ");
      num_added++;
    }

    local_q.pop_front();
  }

  cout << "Added round of " << num_added << " safe Gomory cuts\n";

 CLEANUP:
  if(rval == 1)
    cerr << "Problem in Cut<safeGMI>::add_cut\n";
  return rval;
}
