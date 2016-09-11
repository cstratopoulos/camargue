#ifndef PSEP_SAFEGMI_H
#define PSEP_SAFEGMI_H

#include<vector>
#include<memory>

#include <safemir/src/cplex_slvr.hpp>
#include <safemir/src/ds_cuts.hpp>

#include "lp.h"
#include "cuts.h"

namespace PSEP {
  template<> class Cut<safeGMI> {
  public:
    Cut<safeGMI>(PSEPlp &_m_lp, std::vector<double> &_m_lp_edges,
		 std::vector<int> &_frac_cstat, std::vector<int> &_frac_rstat,
		 std::vector<int> &_support_indices) :
    m_lp(_m_lp), m_lp_edges(_m_lp_edges),
      frac_colstat(_frac_cstat), frac_rowstat(_frac_rstat),
      support_indices(_support_indices) {}

    int cutcall();

  private:
    int init_constraint_info();
    int get_cuts();
    
    PSEPlp &m_lp;
    std::vector<double> &m_lp_edges;
    std::vector<int> &frac_colstat;
    std::vector<int> &frac_rowstat;
    std::vector<int> &support_indices;

    struct SafeMIRGroup {
    SafeMIRGroup(PSEPlp &m_lp):
      vartype(PSEPlp_numcols(&m_lp), 'B'),
	basis_info((SLVRbasisInfo_t *) NULL),
	constraint_matrix((CUTSsystem_t<double> *) NULL),
	tab_row_sparse((CUTSsprow_t<double> *) NULL),
	var_info((CUTSvarInfo_t<double> *) NULL),
	flips((CUTSflips_t *) NULL){
	lp_obj.env = m_lp.cplex_env;
	lp_obj.prob = m_lp.cplex_lp;
	lp_obj.ctype = &vartype[0];
      }
      
      ~SafeMIRGroup(){
	SLVRfreeBasisInfo(&basis_info);
	CUTSfreeSystem(&constraint_matrix);
	CUTSfreeSpRow(&tab_row_sparse);
	CUTSfreeVarInfo(&var_info);
	if(flips->direction) free(flips->direction); free(flips);
      }
      
      std::vector<char> vartype;
      
      SLVRcplex_t lp_obj;
      SLVRbasisInfo_t *basis_info;

      CUTSsystem_t<double> *constraint_matrix;
      
      CUTSsprow_t<double> *tab_row_sparse;
      CUTSsprow_t<double> *current_cut_sparse;
      CUTSsprow_t<double> *best_cut_sparse;

      CUTSvarInfo_t<double> *var_info;
      CUTSflips_t *flips;
    };

    std::unique_ptr<SafeMIRGroup> safe_mir_data;
  };
}

#endif
