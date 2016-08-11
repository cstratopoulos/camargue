#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

#include<utility>
#include<unordered_map>

namespace PSEP {
  enum class SolutionProtocol {
    PURECUT, ABC
      };
  
  namespace LP {
    constexpr double EPSILON = 0.000001;
    constexpr long long DEFAULT_ITLIM = 9223372036800000000;

    enum class Pricing {
      Devex, SlackSteepest, Steepest
    };

    enum class PivType {
      Frac, Subtour, Tour, FathomedTour
    };

    struct Prefs {
    Prefs() : price_method(Pricing::Devex), dp_threshold(-1) {}
    Prefs(Pricing _price, int _dp_threshold) :
      price_method(_price), dp_threshold(_dp_threshold) {}
      
      Pricing price_method;
      int dp_threshold;
    };
  }
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
