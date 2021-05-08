#ifndef _KB_SOLVERS_STARTERS_H_
#define _KB_SOLVERS_STARTERS_H_

#include <stdlib.h>

#include "kepler_benchmarks/helpers.h"

inline double kb_starter_simple(const double M, const double e) {
  return M + sign(sin(M)) * 0.85 * e;
}

inline double kb_starter_markley(const double M, const double e) {
  // M must be in the range [0, pi)
  const double ome = 1. - e;
  const double FACTOR1 = 3 * M_PI / (M_PI - 6 / M_PI);
  const double FACTOR2 = 1.6 / (M_PI - 6 / M_PI);
  const double M2 = M * M;
  const double alpha = FACTOR1 + FACTOR2 * (M_PI - M) / (1 + e);
  const double d = 3 * ome + alpha * e;
  const double alphad = alpha * d;
  const double r = (3 * alphad * (d - ome) + M2) * M;
  const double q = 2 * alphad * ome - M2;
  const double q2 = q * q;
  double w = cbrt(fabs(r) + sqrt(q2 * q + r * r));
  w *= w;
  return (2 * r * w / (w * w + w * q + q2) + M) / d;
}

// Use the second-order series expanion in Raposo-Pulido & Pelaez
// (2017) in the singular corner (eccentricity close to 1, mean
// anomaly close to zero).
inline double kb_starter_rpp(const double M, const double e) {
  const double ome = 1. - e;
  const double sqrt_ome = sqrt(ome);

  const double chi = M / (sqrt_ome * ome);
  const double Lam = sqrt(8 + 9 * chi * chi);
  const double S = cbrt(Lam + 3 * chi);
  const double sigma = 6 * chi / (2 + S * S + 4. / (S * S));
  const double s2 = sigma * sigma;
  const double s4 = s2 * s2;

  const double denom = 1.0 / (s2 + 2);
  const double E =
      sigma * (1 + s2 * ome * denom *
                       ((s2 + 20) / 60. +
                        s2 * ome * denom * denom * (s2 * s4 + 25 * s4 + 340 * s2 + 840) / 1400));

  return E * sqrt_ome;
}

#endif
