#include "util.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

namespace CMR {

LP::Prefs::Prefs() :
  price_method(LP::Pricing::Devex),
  dp_threshold(-1),
  max_per_round(4),
  q_max_size(150) {}

LP::Prefs::Prefs(LP::Pricing _price, int _dp_threshold, int max_round,
		 int q_max) :
  price_method(_price),
  dp_threshold(_dp_threshold),
  max_per_round(max_round),
  q_max_size(q_max) {}


/* zeit function for recording times */

double zeit (void)
{
    struct rusage ru;

    getrusage (RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec) +
           ((double) ru.ru_utime.tv_usec)/1000000.0;
}

double real_zeit (void)
{
    return (double) time (0);
}

std::string LP::piv_string(LP::PivType piv)
{
    switch (piv) {
    case PivType::FathomedTour:
        return "Fathomed tour";
    case PivType::Frac:
        return "Fractional";
    case PivType::Subtour:
        return "Integral subtour";
    case PivType::Tour:
        return "Tour";
    default:
        return "";
    }

    return "";
}

}
