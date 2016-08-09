#ifndef PSEP_GENCUTS_H
#define PSEP_GENCUTS_H

#include<vector>

#include "lp.h"
#include "cuts.h"

namespace PSEP{
  template<> class Cut<general> {
  public:
    Cut<general>(const bool _gomory, const bool _disj, const bool _mir,
		 std::vector<int> &_best_tour_edges, PSEPlp &_m_lp) :
    gencuts(_gomory, _disj, _mir), best_tour_edges(_best_tour_edges),
      m_lp(_m_lp) {}

    int separate();
    int separate(const int edge);

    //private:
    int init_mip();
    int revert_lp();

    int make_all_binary();
    int make_binary(const int edge);

    int check_cuts();
    int num_added(int &frac, int &disj, int &mir);
    
    PSEP::general gencuts;
    std::vector<int> &best_tour_edges;
    PSEPlp &m_lp;
  }
}

#endif
