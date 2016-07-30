#ifndef PSEP_CUTS_H
#define PSEP_CUTS_H

#include <memory>

namespace PSEP {
  template<typename cut_t>
    class Cut {
  public:
    int cut_call(){ return 1;}

  private:
    int separate(){ return 1;}
    int parse_coeffs(){ return 1;}
    int add_cut(){return 1;}

    std::unique_ptr<cut_t> best;
  };
}

#endif
