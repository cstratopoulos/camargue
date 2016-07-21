#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

#include<utility>
#include<unordered_map>

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
PSEP_LP_Prefs() : pricing_choice(0), switching_choice(0),
    dp_threshold(15){}
PSEP_LP_Prefs(int _price, int _switch, int _dp) : pricing_choice(_price),
    switching_choice(_switch), dp_threshold(_dp) {}
  int pricing_choice;
  int switching_choice;
  int dp_threshold;
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

//hash function taken from boost hash_combine



typedef std::pair<int, int> IntPair;
typedef std::unordered_map<IntPair, int> IntPairMap;

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}

double PSEP_zeit (void);
double PSEP_real_zeit (void);

int PSEP_build_xy (int ncount, double *xlist, double *ylist, int gridsize);

bool is_almost_integral(double x);

#endif
