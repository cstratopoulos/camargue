#ifndef PSEP_SAFEGMI_H
#define PSEP_SAFEGMI_H

#include<vector>
#include<memory>

#include<safemir/src/cplex_slvr.hpp>
#include<safemir/src/cutmaster_slvr.hpp>
#include<safemir/src/ds_cuts.hpp>

#include "lp.hpp"
#include "cuts.hpp"

namespace PSEP {
  template<> class Cut<safeGMI> {
  public:
    Cut<safeGMI>(std::vector<int> &_best_tour_edges,
		 PSEPlp &_m_lp, std::vector<double> &_m_lp_edges,
		 std::vector<int> &_frac_cstat,
		 //std::vector<int> &_frac_rstat,
		 std::vector<int> &_support_indices) :
    best_tour_edges(_best_tour_edges),
    m_lp(_m_lp), m_lp_edges(_m_lp_edges),
      frac_colstat(_frac_cstat),
    //frac_rowstat(_frac_rstat),
      support_indices(_support_indices) {}

    int cutcall();

  protected:
    int separate();
    int add_cut();

  private:
    int init_constraint_info();
    int get_tab_rows();

    std::vector<int> &best_tour_edges;
    
    PSEPlp &m_lp;
    std::vector<double> &m_lp_edges;
    std::vector<int> &frac_colstat;
    //std::vector<int> &frac_rowstat;
    std::vector<int> &support_indices;

    struct SafeMIRGroup {
    SafeMIRGroup(PSEPlp &m_lp):
      vartype(PSEPlp_numcols(&m_lp), 'B'),
	basis_info((SLVRbasisInfo_t *) NULL),
	constraint_matrix((CUTSsystem_t<double> *) NULL),
	tableau_rows((CUTSsystem_t<double> *) NULL),
	var_info((CUTSvarInfo_t<double> *) NULL),
	full_x((double *) NULL),
	generated_cuts((CUTSsprowlist_t<double> *) NULL){
	settings.do_mir = true;
	settings.use_log = false;

	settings.max_rank = 1;
	settings.maxcutsperround = 500;

	settings.mir_cutsperfather = 10;
	settings.mir_cutsperround = 1000;
	settings.mir_maxscalings = 15;
	settings.mir_tmin = 1;
	settings.mir_tmax = 3;

	settings.save_cuts = 0;
	
	lp_obj.env = m_lp.cplex_env;
	lp_obj.prob = m_lp.cplex_lp;
	lp_obj.ctype = &vartype[0];
      }
      
      ~SafeMIRGroup(){
	SLVRfreeBasisInfo(&basis_info);
	CUTSfreeSystem(&constraint_matrix);
	CUTSfreeSystem(&tableau_rows);
	CUTSfreeVarInfo(&var_info);
	if(full_x) free(full_x);
	CUTSfreeRowList(&generated_cuts);
      }
      
      std::vector<char> vartype;
      std::vector<double> lp_vranking;

      SLVRcutterSettings_t settings;
      
      SLVRcplex_t lp_obj;
      
      SLVRbasisInfo_t *basis_info;

      CUTSsystem_t<double> *constraint_matrix;
      CUTSsystem_t<double> *tableau_rows;

      CUTSvarInfo_t<double> *var_info;

      double *full_x;

      CUTSsprowlist_t<double> *generated_cuts;
    };

    std::unique_ptr<SafeMIRGroup> safe_mir_data;
  };
}

#endif
