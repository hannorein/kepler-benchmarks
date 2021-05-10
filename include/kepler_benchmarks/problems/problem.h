/**
 * @file  kepler_benchmarks/problems/problem.h
 * @brief The infrastructure for defining and running benchmark problems
 */

#ifndef _KB_PROBLEMS_PROBLEM_H_
#define _KB_PROBLEMS_PROBLEM_H_

#include <stdio.h>
#include <stdlib.h>

#include "kepler_benchmarks/kb_math.h"

/********************************************************************************
 * @fn get_time
 * @brief A cross platform high-resolution timer
 *
 * https://stackoverflow.com/questions/2349776/how-can-i-benchmark-c-code-easily
 *******************************************************************************/
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

/*******************************************************************************
 * @brief The specification of a benchmark problem
 *
 * A "problem" is specified by a single eccentricity and an array of mean
 * anomalies. For each mean anomaly, there is an expected eccentric anomaly and
 * allocated memory for the results of solving Kepler's equation.
 *******************************************************************************/
typedef struct kb_problem {
  long size;         /** The number of mean anomalies in the problem */
  double e;          /** The eccentricity */
  double* M;         /** An array of mean anomalies */
  double* E_expect;  /** The true expected eccentric anomaly */
  double* E_calc;    /** The calculated eccentric anomaly (filled in by the solver) */
  double* sinE_calc; /** The sine of E_calc */
  double* cosE_calc; /** The cosine of E_calc */
} kb_problem;

/***************************************************************************
 * @brief A summary of the results of benchmark run
 *
 * This captures the average runtime and error required to run a particular
 * benchmark.
 **************************************************************************/
typedef struct kb_summary {
  long size;
  double e;
  double runtime;
  double error_max;
  double error_mean;
  double sin_error_max;
  double sin_error_mean;
  double cos_error_max;
  double cos_error_mean;
} kb_summary;

/******************************************************
 * @brief Summarize the results of a benchmark problem
 *
 * @param problem The problem specification/workspace
 * @param runtime The total runtime of the benchmark
 *****************************************************/
kb_summary kb_summarize(const kb_problem* problem, double runtime) {
  kb_summary summary = {0};
  summary.e = problem->e;
  summary.runtime = runtime;
  for (long n = 0; n < problem->size; ++n) {
    double error = problem->E_expect[n] - problem->E_calc[n];
    error = fabs(atan2(sin(error), cos(error)));  // Handle wrapping properly
    double sin_error = fabs(sin(problem->E_expect[n]) - problem->sinE_calc[n]);
    double cos_error = fabs(cos(problem->E_expect[n]) - problem->cosE_calc[n]);
    summary.error_max = max(summary.error_max, error);
    summary.error_mean += error;
    summary.sin_error_max = max(summary.sin_error_max, sin_error);
    summary.sin_error_mean += sin_error;
    summary.cos_error_max = max(summary.cos_error_max, cos_error);
    summary.cos_error_mean += cos_error;
  }
  summary.error_mean /= problem->size;
  summary.sin_error_mean /= problem->size;
  summary.cos_error_mean /= problem->size;
  return summary;
}

/*********************************************************************
 * @brief Allocate the memory required for a particular problem setup
 *
 * @param size The size of the mean anomaly array
 ********************************************************************/
kb_problem* kb_problem_alloc(long size) {
  kb_problem* problem = (kb_problem*)malloc(sizeof(kb_problem));
  problem->size = size;
  problem->e = 0.;
  problem->M = (double*)malloc(size * sizeof(double));
  problem->E_expect = (double*)malloc(size * sizeof(double));
  problem->E_calc = (double*)malloc(size * sizeof(double));
  problem->sinE_calc = (double*)malloc(size * sizeof(double));
  problem->cosE_calc = (double*)malloc(size * sizeof(double));
  return problem;
}

/*****************************************************
 * @brief Free a previously allocated problem setup
 *
 * @param problem The problem specification/workspace
 ****************************************************/
void kb_problem_free(kb_problem* problem) {
  free(problem->M);
  free(problem->E_expect);
  free(problem->E_calc);
  free(problem->sinE_calc);
  free(problem->cosE_calc);
  free(problem);
}

/******************************************************************************
 * @brief Run a benchmark problem and return the summary
 *
 * @param problem      The problem specification/workspace
 * @param solver       The function for solving Kepler's equation
 * @param solver_alloc The (possibly `NULL`) allocation/setup function for the
 *                     solver
 * @param solver_free  The (possibly `NULL`) cleanup function for the solver
 *****************************************************************************/
kb_summary kb_problem_run(kb_problem* problem,
                          double (*solver)(const double, const double, const void* const, double*,
                                           double*),
                          void* (*solver_alloc)(const double), void (*solver_free)(void*)) {
  double e = problem->e;
  double start = 0., end = 0.;
  void* opaque = NULL;

  start = get_time();
  for (long n = 0; n < problem->size; ++n) {
    if (solver_alloc) opaque = solver_alloc(e);
    problem->E_calc[n] =
        solver(problem->M[n], e, opaque, problem->sinE_calc + n, problem->cosE_calc + n);
    if (solver_free) solver_free(opaque);
  }
  end = get_time();

  return kb_summarize(problem, end - start);
}

/********************************************************************************
 * @brief Run a benchmark problem at fixed eccentricity and return the summary
 *
 * This is identical to the `kb_problem_run` function, except it factors out the
 * calls to `solver_alloc` and `solver_free` which some solvers can take
 * advantage of to frontload some calculations at fixed eccentricity.
 *
 * @param problem      The problem specification/workspace
 * @param solver       The function for solving Kepler's equation
 * @param solver_alloc The (possibly `NULL`) allocation/setup function for the
 *                     solver
 * @param solver_free  The (possibly `NULL`) cleanup function for the solver
 *******************************************************************************/
kb_summary kb_problem_run_pre(kb_problem* problem,
                              double (*solver)(const double, const double, const void* const,
                                               double*, double*),
                              void* (*solver_alloc)(const double), void (*solver_free)(void*)) {
  double e = problem->e;
  double start = 0., end = 0.;
  void* opaque = NULL;

  if (solver_alloc) {
    opaque = solver_alloc(e);
  }

  start = get_time();
  for (long n = 0; n < problem->size; ++n) {
    problem->E_calc[n] =
        solver(problem->M[n], e, opaque, problem->sinE_calc + n, problem->cosE_calc + n);
  }
  end = get_time();

  if (solver_free) {
    solver_free(opaque);
  }

  return kb_summarize(problem, end - start);
}

/********************************************************************************
 * @brief Run a benchmark problem over a grid of eccentricities
 *
 * This function will call out to `kb_problem_run` or `kb_problem_run_pre`
 * (depending on the `pre` flag) to actually run the benchmark at a particular
 * eccentricity. The results will be printed to stdout.
 *
 * @param pre          A flag indicating if the pre-computation version of the
 *                     solver should be used
 * @param size         The number of anomalies to run per eccentricity
 * @param e_step       The step size in eccentricity to use when construcint the
 *                     grid
 * @param setup        The problem setup function
 * @param solver       The function for solving Kepler's equation
 * @param solver_alloc The (possibly `NULL`) allocation/setup function for the
 *                     solver
 * @param solver_free  The (possibly `NULL`) cleanup function for the solver
 *******************************************************************************/
void kb_problem_run_ecc_grid(const int pre, const long size, const double e_step,
                             void (*setup)(kb_problem*),
                             double (*solver)(const double, const double, const void* const,
                                              double*, double*),
                             void* (*solver_alloc)(const double), void (*solver_free)(void*)) {
  printf("size,e,runtime,error_max,error_mean,sin_error_max,sin_error_mean,");
  printf("cos_error_max,cos_error_mean\n");
  kb_summary summary;
  kb_problem* problem = kb_problem_alloc(size);
  for (double e = 0.; e < 1; e += e_step) {
    problem->e = e;
    setup(problem);
    if (pre) {
      summary = kb_problem_run(problem, solver, solver_alloc, solver_free);
    } else {
      summary = kb_problem_run_pre(problem, solver, solver_alloc, solver_free);
    }
    printf("%ld,", summary.size);
    printf("%.10e,", summary.e);
    printf("%.6e,", summary.runtime);
    printf("%.6e,", summary.error_max);
    printf("%.6e,", summary.error_mean);
    printf("%.6e,", summary.sin_error_max);
    printf("%.6e,", summary.sin_error_mean);
    printf("%.6e,", summary.cos_error_max);
    printf("%.6e\n", summary.cos_error_mean);
  }
  kb_problem_free(problem);
}

#endif
