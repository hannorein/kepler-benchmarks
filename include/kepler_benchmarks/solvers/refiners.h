#ifndef _KB_SOLVERS_REFINERS_H_
#define _KB_SOLVERS_REFINERS_H_

#include <math.h>

#include "kepler_benchmarks/helpers.h"

double kb_refine_nijenhuis(const double M, const double e, const double E, double* sinE,
                           double* cosE) {
  const double ome = 1. - e;
  const double sE = E - sin(E);
  const double cE = 1. - cos(E);

  const double f_0 = e * sE + E * ome - M;
  const double f_1 = e * cE + ome;
  const double f_2 = e * (E - sE);
  const double f_3 = 1 - f_1;
  const double d_3 = -f_0 / (f_1 - 0.5 * f_0 * f_2 / f_1);
  const double d_4 = -f_0 / (f_1 + 0.5 * d_3 * f_2 + (d_3 * d_3) * f_3 / 6);
  const double d_42 = d_4 * d_4;
  const double dE = -f_0 / (f_1 + 0.5 * d_4 * f_2 + d_4 * d_4 * f_3 / 6 - d_42 * d_4 * f_2 / 24);

  const double EA = E + dE;
  *sinE = sin(EA);
  *cosE = cos(EA);
  return EA;
}

double kb_refine_brandt(const double M, const double e, const double E, double* sinE,
                        double* cosE) {
  double num, denom, dEA, inv_e = 1. / e, sE, cE;

  // Sine and cosine initial guesses using series
  if (E < 0.25 * M_PI) {
    sE = shortsin(E);
    cE = sqrt(1 - sE * sE);
  } else if (E > 0.75 * M_PI) {
    sE = shortsin(M_PI - E);
    cE = -sqrt(1 - sE * sE);
  } else {
    cE = shortsin(0.5 * M_PI - E);
    sE = sqrt(1 - cE * cE);
  }

  // Halley's method to update E
  num = (M - E) * inv_e + sE;
  denom = inv_e - cE;
  dEA = num * denom / (denom * denom + 0.5 * sE * num);

  // Use series to update sin and cos
  if (e < 0.78 || M > 0.4) {
    const double factor = 1 - 0.5 * dEA * dEA;
    *sinE = sE * factor + dEA * cE;
    *cosE = cE * factor - dEA * sE;
  } else {
    // Use Householder's third order method to guarantee performance in the
    // singular corners
    const double factor1 = 1 - 0.5 * dEA * dEA, factor2 = dEA * (1 - dEA * dEA / 6.);
    dEA = num / (denom + dEA * (0.5 * sE + cE * dEA / 6.));
    *sinE = sE * factor1 + cE * factor2;
    *cosE = cE * factor1 - sE * factor2;
  }

  return E + dEA;
}

#endif
