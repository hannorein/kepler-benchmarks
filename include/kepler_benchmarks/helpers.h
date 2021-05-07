#ifndef _KB_HELPERS_H_
#define _KB_HELPERS_H_

#include <math.h>

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

inline int sign(double x) { return (x > 0) - (x < 0); }

double mod_2pi(const double x) {
  const double twopi = 2 * M_PI;
  if (x < twopi && x >= 0) return x;

  if (x > twopi) {
    double M = x - twopi;
    if (M > twopi)
      return fmod(M, twopi);
    else
      return M;
  } else {
    double M = x + twopi;
    if (M < 0)
      return fmod(M, twopi) + twopi;
    else
      return M;
  }
}

#endif
