#ifndef _KB_MATH_H_
#define _KB_MATH_H_

#define _USE_MATH_DEFINES
#include <math.h>

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

inline int sign(double x) { return (x > 0) - (x < 0); }

// Evaluate sine with a series expansion.  We can guarantee that the
// argument will be <=pi/4, and this reaches double precision (within
// a few machine epsilon) at a significantly lower cost than the
// function call to sine that obeys the IEEE standard.
inline double shortsin(const double x) {
  const double if3 = 1. / 6;
  const double if5 = 1. / (6. * 20);
  const double if7 = 1. / (6. * 20 * 42);
  const double if9 = 1. / (6. * 20 * 42 * 72);
  const double if11 = 1. / (6. * 20 * 42 * 72 * 110);
  const double if13 = 1. / (6. * 20 * 42 * 72 * 110 * 156);
  const double if15 = 1. / (6. * 20 * 42 * 72 * 110 * 156 * 210);
  const double x2 = x * x;
  return x *
         (1 - x2 * (if3 -
                    x2 * (if5 - x2 * (if7 - x2 * (if9 - x2 * (if11 - x2 * (if13 - x2 * if15)))))));
}

double mod_2pi(const double x) {
  const double twopi = 2 * M_PI;
  if (x < twopi && x >= 0) return x;

  if (x >= twopi) {
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
