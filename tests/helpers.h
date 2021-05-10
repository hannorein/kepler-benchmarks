#ifndef _KB_TESTS_HELPERS_H_
#define _KB_TESTS_HELPERS_H_

#include <stdio.h>

#include "kepler_benchmarks/kb_math.h"

int run_single_test(const double E, const double e, const double tol,
                    double (*solver)(const double, const double, const void* const, double*,
                                     double*),
                    void* (*solver_alloc)(const double), void (*solver_free)(void*)) {
  double sE = 0., cE = 0.;
  double M = E - e * sin(E);
  void* opaque = NULL;
  if (solver_alloc) opaque = solver_alloc(e);
  double calc = solver(M, e, opaque, &sE, &cE);
  if (solver_free) solver_free(opaque);

  double resid = calc - E;
  resid = fabs(atan2(sin(resid), cos(resid)));  // Handle wrapping properly
  if (resid > tol) {
    printf("Solve failed for E = %.6f, e = %.6f; with error = %.6e\n", E, e, resid / M_PI);
    return 1;
  }

  resid = fabs(sE - sin(E));
  if (resid > tol) {
    printf("Solve failed for sin(E) = %.6f, e = %.6f; with error = %.6e\n", sin(E), e, resid);
    return 1;
  }

  resid = fabs(cE - cos(E));
  if (resid > tol) {
    printf("Solve failed for cos(E) = %.6f, e = %.6f; with error = %.6e\n", cos(E), e, resid);
    return 1;
  }

  return 0;
}

int run_solver_test(const double tol,
                    double (*solver)(const double, const double, const void* const, double*,
                                     double*),
                    void* (*solver_alloc)(const double), void (*solver_free)(void*)) {
  // Test high eccentricity at some sore spots
  if (run_single_test(0., 1. - 1e-6, tol, solver, solver_alloc, solver_free)) return 1;
  if (run_single_test(2 * M_PI, 1. - 1e-6, tol, solver, solver_alloc, solver_free)) return 1;
  if (run_single_test(-226.2, 1. - 1e-6, tol, solver, solver_alloc, solver_free)) return 1;
  if (run_single_test(-170.4, 0.9939879759519037, tol, solver, solver_alloc, solver_free))
    return 1;

  // Run some tests where E is exactly pi, 2*pi, -pi, and -2*pi
  for (double e = 0.0; e < 1.0; e += 0.0134) {
    if (run_single_test(M_PI, e, tol, solver, solver_alloc, solver_free)) return 1;
    if (run_single_test(-M_PI, e, tol, solver, solver_alloc, solver_free)) return 1;
    if (run_single_test(2 * M_PI, e, tol, solver, solver_alloc, solver_free)) return 1;
    if (run_single_test(-2 * M_PI, e, tol, solver, solver_alloc, solver_free)) return 1;
  }

  // Then finally a small grid
  for (double e = 0.0; e < 1.0; e += 0.0134) {
    for (double E = -10.2; E < 12.0; E += 0.178) {
      if (run_single_test(E, e, tol, solver, solver_alloc, solver_free)) return 1;
    }
  }

  return 0;
}

#endif
