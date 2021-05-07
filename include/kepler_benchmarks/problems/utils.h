#ifndef _KB_PROBLEMS_UTILS_H_
#define _KB_PROBLEMS_UTILS_H_

#include <math.h>
#include <stdlib.h>

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

// Cross platform high-resolution timer
// Ref: https://stackoverflow.com/questions/2349776/how-can-i-benchmark-c-code-easily
#ifdef WIN32

#include <windows.h>
double get_time() {
  LARGE_INTEGER t, f;
  QueryPerformanceCounter(&t);
  QueryPerformanceFrequency(&f);
  return (double)t.QuadPart / (double)f.QuadPart;
}

#else

#include <sys/resource.h>
#include <sys/time.h>

double get_time() {
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec * 1e-6;
}

#endif

typedef struct kb_benchmark_result {
  double runtime;
  double error_max;
  double error_mean;
  double sin_error_max;
  double sin_error_mean;
  double cos_error_max;
  double cos_error_mean;
} kb_test_result;

typedef struct kb_benchmark_problem {
  int fixed_ecc;
  long size;
  double* e;
  double* M;
  double* E_expect;
  double* E_calc;
  double* sinE_calc;
  double* cosE_calc;
} kb_test;

kb_test* kb_test_alloc(long size) {
  kb_test* test = (kb_test*)malloc(sizeof(kb_test));
  test->fixed_ecc = 0;
  test->size = size;
  test->e = (double*)malloc(size * sizeof(double));
  test->M = (double*)malloc(size * sizeof(double));
  test->E_expect = (double*)malloc(size * sizeof(double));
  test->E_calc = (double*)malloc(size * sizeof(double));
  test->sinE_calc = (double*)malloc(size * sizeof(double));
  test->cosE_calc = (double*)malloc(size * sizeof(double));
  return test;
}

void kb_test_free(kb_test* test) {
  free(test->e);
  free(test->M);
  free(test->E_expect);
  free(test->E_calc);
  free(test->sinE_calc);
  free(test->cosE_calc);
  free(test);
}

kb_test_result kb_test_run(kb_test* test,
                           double (*solver)(const double, const double, const void* const, double*,
                                            double*),
                           void* (*solver_alloc)(const double), void (*solver_free)(void*)) {
  double start = 0., end = 0.;
  void* precalc = NULL;  // precalc can be nothing

  if (test->fixed_ecc) {
    double e = test->e[0];
    if (solver_alloc) {
      precalc = solver_alloc(e);
    }

    start = get_time();
    for (long n = 0; n < test->size; ++n) {
      test->E_calc[n] = solver(test->M[n], e, precalc, test->sinE_calc + n, test->cosE_calc + n);
    }
    end = get_time();

    if (solver_free) {
      solver_free(precalc);
    }
  } else {
    double e;
    start = get_time();
    for (long n = 0; n < test->size; ++n) {
      e = test->e[n];
      if (solver_alloc) precalc = solver_alloc(e);
      test->E_calc[n] = solver(test->M[n], e, precalc, test->sinE_calc + n, test->cosE_calc + n);
      if (solver_free) solver_free(precalc);
    }
    end = get_time();
  }

  // Compute the results
  kb_test_result result = {0};
  result.runtime = end - start;
  for (long n = 0; n < test->size; ++n) {
    double error = fabs(test->E_expect[n] - test->E_calc[n]);
    double sin_error = fabs(sin(test->E_expect[n]) - test->sinE_calc[n]);
    double cos_error = fabs(cos(test->E_expect[n]) - test->cosE_calc[n]);
    result.error_max = max(result.error_max, error);
    result.error_mean += error;
    result.sin_error_max = max(result.sin_error_max, sin_error);
    result.sin_error_mean += sin_error;
    result.cos_error_max = max(result.cos_error_max, cos_error);
    result.cos_error_mean += cos_error;
  }
  result.error_mean /= test->size;
  result.sin_error_mean /= test->size;
  result.cos_error_mean /= test->size;
  return result;
}

#endif
