#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

namespace LP {
  const double EPSILON = 0.000001;
  const long DEFAULT_ITLIM = 9223372036800000000;

  namespace PRICING {
    const int DEVEX = 0;
    const int STEEPEST = 1;
    const int STEEPEST_REAL = 2;

    static int choice;

    namespace SWITCHING{
      const int OFF = 0;
      const int DYNAMIC = 1;
      const int START = 2;

      static int choice;
    }
  }
}

namespace UTIL {
  static int seed = 0;
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
