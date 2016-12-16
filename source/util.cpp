#include "util.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

namespace CMR {


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
