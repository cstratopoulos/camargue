#ifndef PSEP_SAFEGMI_H
#define PSEP_SAFEGMI_H

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

    int test();

  private:
    PSEPlp &m_lp;
    std::vector<double> &m_lp_edges;
    std::vector<int> &frac_colstat;
    std::vector<int> &frac_rowstat;
    std::vector<int> &support_indices;
  };
}

#endif
