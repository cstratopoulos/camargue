/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief PRIMAL SEPARATION OF SAFE GOMORY MIXED-INTEGER CUTS
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_SAFEGMI_H
#define CMR_SAFEGMI_H

#include<vector>
#include<memory>

#include<safemir/src/cplex_slvr.hpp>
#include<safemir/src/cutmaster_slvr.hpp>
#include<safemir/src/ds_cuts.hpp>

#include "lp.hpp"
#include "cuts.hpp"

namespace CMR {

/** Wrapper class for using cuts from the safemir code. */
struct cut_obj {
  cut_obj(bool _exact, int _zeros, double _viol,
	  CUTSsprow_t<double> *_row):
    exact(_exact),
    zeros(_zeros),
    viol(_viol),
    row(_row){}

  bool exact; /**< Is the cut tight by double `==`.*/
  int zeros; /**< Number of zeros in the cut. */
  double viol; /** Amount by which cut is violated. */

  CUTSsprow_t<double> *row; /**< Inherited structure from safemir. */
};

/** Idiomatic implementation of lexicographic order on cut_obj.
 * Exactness is most important, then density, then violation.
 */
inline static bool operator >(const cut_obj &a, const cut_obj &b) {
  return std::tie(a.exact, a.zeros, a.viol) >
    std::tie(b.exact, b.zeros, b.viol);
}

/** Class for primal separation of safe GMI cuts.
 * This class performs exact primal separation of safe GMI cuts using the
 * algorithm and numerical safety framework of Cook, Dash, Fukasawa, and 
 * Goycoolea. All the members and classes/structs are wrappers to the code
 * provided online as a supplement to the paper
 *  Numerically safe Gomory mixed-integer
 * cuts, W. Cook,  S. Dash, R. Fukasawa, and M. Goycoolea, 
 * INFORMS Journal on Computing 21 (2009) 641--649.
 */
template<> class Cut<safeGMI> {
public:
  Cut<safeGMI>(std::vector<int> &_best_tour_edges,
	       CMRlp &_m_lp, std::vector<double> &_m_lp_edges,
	       std::vector<int> &_frac_cstat,
	       std::vector<int> &_frac_rstat,
	       std::vector<int> &_support_indices,
	       const int max_per_round) :
    best_tour_edges(_best_tour_edges),
    m_lp(_m_lp), m_lp_edges(_m_lp_edges),
    frac_colstat(_frac_cstat),
    frac_rowstat(_frac_rstat),
    support_indices(_support_indices),
    local_q(max_per_round){}

  int cutcall();

protected:
  int separate();
  int add_cut();

private:
  int init_constraint_info();
  int get_tab_rows();

  std::vector<int> &best_tour_edges;
    
  CMRlp &m_lp;
  std::vector<double> &m_lp_edges;
  std::vector<int> &frac_colstat;
  std::vector<int> &frac_rowstat;
  std::vector<int> &support_indices;

  struct SafeMIRGroup {
    SafeMIRGroup(CMRlp &m_lp):
      vartype(CMRlp_numcols(&m_lp), 'B'),
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

  CutQueue<cut_obj> local_q;
};
}

#endif
