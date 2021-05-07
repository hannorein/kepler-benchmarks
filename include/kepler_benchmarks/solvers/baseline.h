/**
 * @file  kepler_benchmarks/solvers/baseline.h
 * @brief A trivial dummy solver
 *
 * This solver only returns the input value and computes it sine and cosine. It
 * will only give the right answers for zero eccentricity, but it is used to
 * measure the runtime of evaluating trig functions.
 */

#ifndef _KB_SOLVERS_BASELINE_H_
#define _KB_SOLVERS_BASELINE_H_

#include <math.h>

double kb_solver_baseline(const double M, const double e, const void *opaque, double *sinE,
                          double *cosE) {
  (void)e;
  (void)opaque;
  *sinE = sin(M);
  *cosE = cos(M);
  return M;
}

#endif
