#ifndef _KB_SOLVERS_SIMPLE_H_
#define _KB_SOLVERS_SIMPLE_H_

#include <math.h>
#include <stdlib.h>

// Dummy Kepler solver
double kb_solver_simple(const double M, const double e, const void *precalc, double *sinE,
                        double *cosE) {
  double precalculated_value = *(double *)precalc;
  double E = M * e * precalculated_value;
  *sinE = sin(E);
  *cosE = cos(E);
  return E;
}

// Precalculate values
// Pointer can point to array, struct, etc
void *kb_solver_simple_alloc(const double e) {
  (void)e;
  double *p = (double *)malloc(sizeof(double));
  *p = sqrt(1.);
  return p;
}

// Free precalculated values
void kb_solver_simple_free(void *p) { free(p); }

#endif
