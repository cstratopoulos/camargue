#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

namespace LP {
  static const double EPSILON = 0.000001;
  static const long DEFAULT_ITLIM = 9223372036800000000;

  namespace PRICING {
    static const int DEVEX = 0;
    static const int STEEPEST = 1;
    static const int STEEPEST_REAL = 2;

    namespace SWITCHING{
      static const int OFF = 0;
      static const int DYNAMIC = 1;
      static const int START = 2;
    }
  }
}

struct PSEP_LP_Prefs {
PSEP_LP_Prefs() : pricing_choice(0), switching_choice(0) {}
PSEP_LP_Prefs(int _price, int _switch) : pricing_choice(_price),
    switching_choice(_switch) {}
  int pricing_choice;
  int switching_choice;
};

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
