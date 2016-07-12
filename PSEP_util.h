#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

namespace LP {
  const double EPSILON = 0.000001;
  const long DEFAULT_ITLIM = 9223372036800000000;
  static bool devex_switch = false;
}

namespace PIVOT {
  const int FRAC = 0;
  const int SUBTOUR = 1;
  const int TOUR = 2;
  const int FATHOMED_TOUR = 3;
}

double PSEP_zeit (void);
double PSEP_real_zeit (void);

int PSEP_build_xy (int ncount, double *xlist, double *ylist, int gridsize);

bool is_almost_integral(double x);

#endif
