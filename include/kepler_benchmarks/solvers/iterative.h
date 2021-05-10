/**
 * @file  kepler_benchmarks/solvers/iterative.h
 * @brief A set of commonly used iterative solvers
 */

#ifndef _KB_SOLVERS_ITERATIVE_H_
#define _KB_SOLVERS_ITERATIVE_H_

#include "kepler_benchmarks/kb_math.h"
#include "kepler_benchmarks/solvers/starters.h"

#ifndef KB_ITER_TOL
#define KB_ITER_TOL (1.234e-12)
#endif

#ifndef KB_MAX_ITER
#define KB_MAX_ITER 30
#endif

/********************************************************************************
 * @brief A first order iterative solver
 *
 * This is based on the implementation in the "batman" code, which is (in turn)
 * based on the algorithm in the algorithm in the Murray & Correia chapter of
 * the "Exoplanets" textbook.
 *******************************************************************************/
double kb_solver_iter_first(const double M, const double e, const void *opaque, double *sinE,
                            double *cosE) {
  (void)opaque;  // unused

  double E = kb_starter_simple(M, e);
  double sE = sin(E), cE = cos(E);
  double fe, fs;

  for (int i = 0; i < KB_MAX_ITER; ++i) {
    fe = E - e * sE - M;
    if (fabs(fe) < KB_ITER_TOL) break;
    fs = 1 - e * cE;
    E = E - fe / fs;
    sE = sin(E);
    cE = cos(E);
  }

  *sinE = sE;
  *cosE = cE;
  return E;
}

/************************************************************
 * @brief A third order iterative solver
 *
 * This is based on the implementation in the "radvel" code.
 ***********************************************************/
double kb_solver_iter_third(const double M, const double e, const void *opaque, double *sinE,
                            double *cosE) {
  (void)opaque;  // unused

  double E = kb_starter_simple(M, e);
  double sE = sin(E), cE = cos(E);
  double fi, fip, fipp, fippp, d1;

  for (int i = 0; i < KB_MAX_ITER; ++i) {
    fi = (E - e * sE - M);
    if (fabs(fi) < KB_ITER_TOL) break;

    fip = 1 - e * cE;
    fipp = e * sE;
    fippp = 1 - fip;

    d1 = -fi / fip;
    d1 = -fi / (fip + d1 * fipp / 2.0);
    d1 = -fi / (fip + d1 * fipp / 2.0 + d1 * d1 * fippp / 6.0);
    E += d1;

    sE = sin(E);
    cE = cos(E);
  }

  *sinE = sE;
  *cosE = cE;
  return E;
}

#endif
