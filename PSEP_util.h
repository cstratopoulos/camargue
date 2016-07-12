#ifndef __PSEP_UTIL_H
#define __PSEP_UTIL_H

#define LP_EPSILON 0.000001
#define LP_DEFAULT_ITLIM 9223372036800000000

double PSEP_zeit (void);
double PSEP_real_zeit (void);

int PSEP_build_xy (int ncount, double *xlist, double *ylist, int gridsize);

bool is_almost_integral(double x);

#endif
