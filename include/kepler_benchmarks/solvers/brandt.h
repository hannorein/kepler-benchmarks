#ifndef _KB_SOLVERS_BRANDT_H_
#define _KB_SOLVERS_BRANDT_H_

#include "kepler_benchmarks/solvers/refiners.h"
#include "kepler_benchmarks/solvers/starters.h"

#ifndef KB_BRANDT_TOL
#define KB_BRANDT_TOL (1.e-12)
#endif
typedef struct kb_solver_brandt_info {
  double e;
  double bounds[13];
  double EA_tab[6 * 13];
} kb_solver_brandt_info;

void* kb_solver_brandt_alloc(const double e) {
  kb_solver_brandt_info* opaque = (kb_solver_brandt_info*)malloc(sizeof(kb_solver_brandt_info));
  opaque->e = e;
  if (e < KB_BRANDT_TOL) return (void*)opaque;

  const double pi = M_PI;
  const double pi_d_12 = M_PI / 12;
  const double pi_d_6 = M_PI / 6;
  const double pi_d_4 = M_PI / 4;
  const double pi_d_3 = M_PI / 3;
  const double fivepi_d_12 = M_PI * 5. / 12;
  const double pi_d_2 = M_PI / 2;
  const double sevenpi_d_12 = M_PI * 7. / 12;
  const double twopi_d_3 = M_PI * 2. / 3;
  const double threepi_d_4 = M_PI * 3. / 4;
  const double fivepi_d_6 = M_PI * 5. / 6;
  const double elevenpi_d_12 = M_PI * 11. / 12;

  const double g2s_e = 0.2588190451025207623489 * e;
  const double g3s_e = 0.5 * e;
  const double g4s_e = 0.7071067811865475244008 * e;
  const double g5s_e = 0.8660254037844386467637 * e;
  const double g6s_e = 0.9659258262890682867497 * e;
  const double g2c_e = g6s_e;
  const double g3c_e = g5s_e;
  const double g4c_e = g4s_e;
  const double g5c_e = g3s_e;
  const double g6c_e = g2s_e;

  opaque->bounds[0] = 0;
  opaque->bounds[1] = pi_d_12 - g2s_e;
  opaque->bounds[2] = pi_d_6 - g3s_e;
  opaque->bounds[3] = pi_d_4 - g4s_e;
  opaque->bounds[4] = pi_d_3 - g5s_e;
  opaque->bounds[5] = fivepi_d_12 - g6s_e;
  opaque->bounds[6] = pi_d_2 - e;
  opaque->bounds[7] = sevenpi_d_12 - g6s_e;
  opaque->bounds[8] = twopi_d_3 - g5s_e;
  opaque->bounds[9] = threepi_d_4 - g4s_e;
  opaque->bounds[10] = fivepi_d_6 - g3s_e;
  opaque->bounds[11] = elevenpi_d_12 - g2s_e;
  opaque->bounds[12] = pi;

  double x;

  opaque->EA_tab[1] = 1. / (1. - e);
  opaque->EA_tab[2] = 0;

  x = 1. / (1 - g2c_e);
  opaque->EA_tab[7] = x;
  opaque->EA_tab[8] = -0.5 * g2s_e * x * x * x;
  x = 1. / (1 - g3c_e);
  opaque->EA_tab[13] = x;
  opaque->EA_tab[14] = -0.5 * g3s_e * x * x * x;
  x = 1. / (1 - g4c_e);
  opaque->EA_tab[19] = x;
  opaque->EA_tab[20] = -0.5 * g4s_e * x * x * x;
  x = 1. / (1 - g5c_e);
  opaque->EA_tab[25] = x;
  opaque->EA_tab[26] = -0.5 * g5s_e * x * x * x;
  x = 1. / (1 - g6c_e);
  opaque->EA_tab[31] = x;
  opaque->EA_tab[32] = -0.5 * g6s_e * x * x * x;

  opaque->EA_tab[37] = 1;
  opaque->EA_tab[38] = -0.5 * e;

  x = 1. / (1 + g6c_e);
  opaque->EA_tab[43] = x;
  opaque->EA_tab[44] = -0.5 * g6s_e * x * x * x;
  x = 1. / (1 + g5c_e);
  opaque->EA_tab[49] = x;
  opaque->EA_tab[50] = -0.5 * g5s_e * x * x * x;
  x = 1. / (1 + g4c_e);
  opaque->EA_tab[55] = x;
  opaque->EA_tab[56] = -0.5 * g4s_e * x * x * x;
  x = 1. / (1 + g3c_e);
  opaque->EA_tab[61] = x;
  opaque->EA_tab[62] = -0.5 * g3s_e * x * x * x;
  x = 1. / (1 + g2c_e);
  opaque->EA_tab[67] = x;
  opaque->EA_tab[68] = -0.5 * g2s_e * x * x * x;

  opaque->EA_tab[73] = 1. / (1. + e);
  opaque->EA_tab[74] = 0;

  double B0, B1, B2, idx;
  int i, k;
  for (i = 0; i < 12; i++) {
    idx = 1. / (opaque->bounds[i + 1] - opaque->bounds[i]);
    k = 6 * i;
    opaque->EA_tab[k] = i * pi_d_12;

    B0 = idx * (-opaque->EA_tab[k + 2] - idx * (opaque->EA_tab[k + 1] - idx * pi_d_12));
    B1 =
        idx * (-2 * opaque->EA_tab[k + 2] - idx * (opaque->EA_tab[k + 1] - opaque->EA_tab[k + 7]));
    B2 = idx * (opaque->EA_tab[k + 8] - opaque->EA_tab[k + 2]);

    opaque->EA_tab[k + 3] = B2 - 4 * B1 + 10 * B0;
    opaque->EA_tab[k + 4] = (-2 * B2 + 7 * B1 - 15 * B0) * idx;
    opaque->EA_tab[k + 5] = (B2 - 3 * B1 + 6 * B0) * idx * idx;
  }

  return (void*)opaque;
}

void kb_solver_brandt_free(void* opaque) { free(opaque); }

double kb_solver_brandt_starter_fixed_e(const double M, const void* opaque) {
  const kb_solver_brandt_info* info = (const kb_solver_brandt_info*)opaque;
  const double e = info->e;
  const double* bounds = info->bounds;
  const double* EA_tab = info->EA_tab;

  int j, k;
  double EA, dx;
  if (e < 0.78 || 2 * M + (1 - e) > 0.2) {
    for (j = 11; j > 0; --j)
      if (M > bounds[j]) break;

    k = 6 * j;
    dx = M - bounds[j];
    EA =
        EA_tab[k] + dx * (EA_tab[k + 1] +
                          dx * (EA_tab[k + 2] +
                                dx * (EA_tab[k + 3] + dx * (EA_tab[k + 4] + dx * EA_tab[k + 5]))));
  } else {
    EA = kb_starter_rpp(M, e);
  }
  return EA;
}

double kb_solver_brandt_starter(const double M, const double e) {
  double bounds[13];
  double EA_tab[9];

  int k;
  double MA = M, EA, x, y;
  double B0, B1, B2, dx, idx;

  // Series expansion
  if (2 * MA + 1 - e < 0.2) {
    EA = kb_starter_rpp(MA, e);
  } else {
    const double pi_d_12 = M_PI / 12;
    const double pi_d_6 = M_PI / 6;
    const double pi_d_4 = M_PI / 4;
    const double pi_d_3 = M_PI / 3;
    const double fivepi_d_12 = M_PI * 5. / 12;
    const double pi_d_2 = M_PI / 2;
    const double sevenpi_d_12 = M_PI * 7. / 12;
    const double twopi_d_3 = M_PI * 2. / 3;
    const double threepi_d_4 = M_PI * 3. / 4;
    const double fivepi_d_6 = M_PI * 5. / 6;
    const double elevenpi_d_12 = M_PI * 11. / 12;

    const double g2s_e = 0.2588190451025207623489 * e;
    const double g3s_e = 0.5 * e;
    const double g4s_e = 0.7071067811865475244008 * e;
    const double g5s_e = 0.8660254037844386467637 * e;
    const double g6s_e = 0.9659258262890682867497 * e;

    // Polynomial boundaries given in Raposo-Pulido & Pelaez
    bounds[0] = 0;
    bounds[1] = pi_d_12 - g2s_e;
    bounds[2] = pi_d_6 - g3s_e;
    bounds[3] = pi_d_4 - g4s_e;
    bounds[4] = pi_d_3 - g5s_e;
    bounds[5] = fivepi_d_12 - g6s_e;
    bounds[6] = pi_d_2 - e;
    bounds[7] = sevenpi_d_12 - g6s_e;
    bounds[8] = twopi_d_3 - g5s_e;
    bounds[9] = threepi_d_4 - g4s_e;
    bounds[10] = fivepi_d_6 - g3s_e;
    bounds[11] = elevenpi_d_12 - g2s_e;
    bounds[12] = M_PI;

    // Which interval?
    for (k = 11; k > 0; k--) {
      if (MA > bounds[k]) break;
    }
    // if (k < 0) k = 0;

    // Values at the two endpoints.
    EA_tab[0] = k * pi_d_12;
    EA_tab[6] = (k + 1) * pi_d_12;

    // First two derivatives at the endpoints. Left endpoint first.
    int sign = (k >= 6) ? 1 : -1;

    x = 1 / (1 - ((6 - k) * pi_d_12 + sign * bounds[abs(6 - k)]));
    y = -0.5 * (k * pi_d_12 - bounds[k]);
    EA_tab[1] = x;
    EA_tab[2] = y * x * x * x;

    x = 1 / (1 - ((5 - k) * pi_d_12 + sign * bounds[abs(5 - k)]));
    y = -0.5 * ((k + 1) * pi_d_12 - bounds[k + 1]);
    EA_tab[7] = x;
    EA_tab[8] = y * x * x * x;

    // Solve a matrix equation to get the rest of the coefficients.
    idx = 1 / (bounds[k + 1] - bounds[k]);

    B0 = idx * (-EA_tab[2] - idx * (EA_tab[1] - idx * pi_d_12));
    B1 = idx * (-2 * EA_tab[2] - idx * (EA_tab[1] - EA_tab[7]));
    B2 = idx * (EA_tab[8] - EA_tab[2]);

    EA_tab[3] = B2 - 4 * B1 + 10 * B0;
    EA_tab[4] = (-2 * B2 + 7 * B1 - 15 * B0) * idx;
    EA_tab[5] = (B2 - 3 * B1 + 6 * B0) * idx * idx;

    // Now use the coefficients of this polynomial to get the initial guess.
    dx = MA - bounds[k];
    EA =
        EA_tab[0] +
        dx * (EA_tab[1] + dx * (EA_tab[2] + dx * (EA_tab[3] + dx * (EA_tab[4] + dx * EA_tab[5]))));
  }
  return EA;
}

double kb_solver_brandt(const double M, const double e, const void* opaque, double* sinE,
                        double* cosE) {
  (void)(opaque);

  if (e < KB_BRANDT_TOL) {
    *sinE = sin(M);
    *cosE = cos(M);
    return M;
  }

  int MAsign = 1;
  double MA = mod_2pi(M);
  printf("M = %f; MA = %f\n", M, MA);
  if (MA > M_PI) {
    MAsign = -1;
    MA = 2 * M_PI - MA;
  }

  double E = kb_solver_brandt_starter(MA, e);
  E = MAsign * kb_refine_brandt(MA, e, E, sinE, cosE);
  *sinE *= MAsign;
  return E;
}

double kb_solver_brandt_fixed_e(const double M, const double e, const void* opaque, double* sinE,
                                double* cosE) {
  if (e < KB_BRANDT_TOL) {
    *sinE = sin(M);
    *cosE = cos(M);
    return M;
  }

  int MAsign = 1;
  double MA = mod_2pi(M);
  if (MA > M_PI) {
    MAsign = -1;
    MA = 2 * M_PI - MA;
  }

  double E = kb_solver_brandt_starter_fixed_e(MA, opaque);
  E = MAsign * kb_refine_brandt(MA, e, E, sinE, cosE);
  *sinE *= MAsign;
  return E;
}

#endif
