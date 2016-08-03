#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

#include<utility>
#include<unordered_map>

namespace LP {
  constexpr double EPSILON = 0.000001;
  constexpr long long DEFAULT_ITLIM = 9223372036800000000;

  namespace PRICING {
    constexpr int DEVEX = 0;
    constexpr int STEEPEST = 1;
    constexpr int STEEPEST_REAL = 2;

    namespace SWITCHING{
      constexpr int OFF = 0;
      constexpr int DYNAMIC = 1;
      constexpr int START = 2;
    }
  }
}

struct PSEP_LP_Prefs {
PSEP_LP_Prefs() : pricing_choice(0), switching_choice(0),
    dp_threshold(15), jumpstart(false),
    redcost_fixing(true){}
PSEP_LP_Prefs(int _price, int _switch, int _dp,
	      bool _jump, bool _redfix, int _max) : pricing_choice(_price),
    switching_choice(_switch), dp_threshold(_dp), jumpstart(_jump),
    redcost_fixing(_redfix){}
  int pricing_choice;
  int switching_choice;
  int dp_threshold;
  bool jumpstart;
  bool redcost_fixing;
};

namespace UTIL {
  static int seed = 0;
}

namespace PSEP{
typedef enum class PivType {
  FRAC, SUBTOUR, TOUR, FATHOMED_TOUR
    } PivType;
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

#endif
