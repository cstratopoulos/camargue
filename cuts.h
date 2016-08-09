#ifndef PSEP_CUTS_H
#define PSEP_CUTS_H

#include <memory>
#include <vector>

#include "tooth.h"

namespace PSEP {
  
  struct seg {
  seg(int _start, int _end, double _slack) :
    start(_start), end(_end), viol(_slack) {}
    int start;
    int end;
    double viol;

    bool operator <(const seg &val) const {
      return viol > val.viol;
    }
  };

  struct blossom {
  blossom(std::vector<int> &_handle, int _cut_edge, double _val) :
    handle(_handle), cut_edge(_cut_edge), cut_val(_val){}

    bool operator< (const blossom &val) const {
      return cut_val < val.cut_val;
    }

    std::vector<int> handle;
    int cut_edge;
    double cut_val;
  };

  struct domino {
    domino(){}
    domino(std::vector<int> &_handle,
	   std::vector<std::shared_ptr<PSEP_CandTooth::SimpleTooth> > _teeth) :
    handle(_handle), used_teeth(_teeth) {}

    std::vector<int> handle;
    std::vector<std::shared_ptr<PSEP_CandTooth::SimpleTooth> > used_teeth;
  };

  struct general {
  general(bool _gomory, bool _disj, bool _mir) :
    added_range(IntPair(0,0)), branch_var(-1), gomory_frac(_gomory),
      disjunctive(_disj), rounding(_mir) {}
  general(const int _var) :
    added_range(IntPair(0,0)), branch_var(_var), gomory_frac(false),
      disjunctive(true), rounding(false) {}
    
    IntPair added_range;
    const int branch_var;

    bool gomory_frac;
    bool disjunctive;
    bool rounding;
  };

  
  template<typename cut_t>
    class Cut {
  /* public: */
  /*   int cut_call(){ return 1;} */

  /* private: */
  /*   int separate(){ return 1;} */
  /*   int parse_coeffs(){ return 1;} */
  /*   int add_cut(){return 1;} */

  /*   std::unique_ptr<cut_t> best; */
  };
}

#endif
