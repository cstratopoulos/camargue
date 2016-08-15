#ifndef PSEP_GENCUTS_H
#define PSEP_GENCUTS_H

#include<array>
#include<vector>

#include "lp.h"
#include "cuts.h"

namespace PSEP{
  template<> class Cut<general> {
  public:
    Cut<general>(const bool _gomory, const bool _disj, const bool _mir,
		 std::vector<int> &_best_tour_edges, PSEPlp &_m_lp,
		 std::vector<double> &_m_lp_edges,
		 std::vector<int> &_frac_cstat, std::vector<int> &_frac_rstat) :
    gencuts(_gomory, _disj, _mir), best_tour_edges(_best_tour_edges),
      m_lp(_m_lp), m_lp_edges(_m_lp_edges),
      frac_colstat(_frac_cstat), frac_rowstat(_frac_rstat) {}

    int separate(const double piv_val);

  private:
    struct generated_cut {
    generated_cut(): found_cut(false), best_nonzeros(100000),
	is_best_exact(false), best_viol(0){}
      generated_cut(const int numcols, const int numrows,
		    const std::vector<int> &_best_tour_edges,
		    const std::vector<double> &_m_lp_edges,
		    const std::vector<int> &_frac_cstat,
		    const std::vector<int> &_frac_rstat):
      found_cut(false), initial_numrows(numrows), best_nonzeros(numcols),
	is_best_exact(false),
      best_tour_edges(&_best_tour_edges[0]), m_lp_edges(&_m_lp_edges[0]),
	frac_colstat(&_frac_cstat[0]), frac_rowstat(&_frac_rstat[0]),
	best_viol(0.0){
	coefficient_buffer.resize(numcols);
	index_buffer.resize(numcols);

	best_coeffs.resize(numcols);
	best_indices.resize(numcols);
      }

      bool found_cut;

      int initial_numrows;

      std::vector<double> coefficient_buffer;
      std::vector<int> index_buffer;

      std::vector<double> best_coeffs;
      std::vector<int> best_indices;
      int best_nonzeros;
      char best_sense;
      double best_rhs;
      bool is_best_exact;
      
      int const *best_tour_edges;
      double const *m_lp_edges;
      int const *frac_colstat;
      int const *frac_rowstat;
      double best_viol;
    };
    
    int init_mip(const double piv_val, generated_cut &callback_arg);
    int revert_lp();

    int make_all_binary();

    int deletion_row;
    
    PSEP::general gencuts;
    std::vector<int> &best_tour_edges;
    PSEPlp &m_lp;
    std::vector<double> &m_lp_edges;
    std::vector<int> &frac_colstat;
    std::vector<int> &frac_rowstat;

    static int CPXPUBLIC solvecallback(CPXCENVptr env, void *cbdata,
				       int wherefrom, void *cbhandle,
				       int *useraction_p);
  };
}

#endif
