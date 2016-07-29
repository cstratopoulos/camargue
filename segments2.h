#ifndef PSEP_SEGCUTS_2_H
#define PSEP_SEGCUTS_2_H

#include "cuts.h"

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
}

#endif
