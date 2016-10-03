#include "PSEP_util.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

using namespace PSEP;

LP::Prefs::Prefs() :
  price_method(LP::Pricing::Devex),
  dp_threshold(-1),
  max_per_round(1),
  q_max_size(50) {}

LP::Prefs::Prefs(LP::Pricing _price, int _dp_threshold, int max_round,
		 int q_max) :
  price_method(_price),
  dp_threshold(_dp_threshold),
  max_per_round(max_round),
  q_max_size(q_max){}


/* zeit function for recording times */

double PSEP::zeit (void)
{
    struct rusage ru;

    getrusage (RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec) +
           ((double) ru.ru_utime.tv_usec)/1000000.0;
}

double PSEP::real_zeit (void)
{
    return (double) time (0);
}
