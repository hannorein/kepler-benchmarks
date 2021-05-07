/**
 * @file  kepler_benchmarks/problems/basic.h
 * @brief The simplest benchmark problem with a uniform grid test points
 */

#ifndef _KB_PROBLEMS_BASIC_H_
#define _KB_PROBLEMS_BASIC_H_

#include "kepler_benchmarks/problems/problem.h"

/******************************************************************************
 * @brief Setup a simple benchmark problem
 *
 * This problem is a uniform grid in eccentric anomaly in the range [0, 2*pi).
 *****************************************************************************/
void kb_problem_setup_basic(kb_problem *problem) {
  double e = problem->e;
  for (long n = 0; n < problem->size; ++n) {
    double E = 2 * M_PI * n / (problem->size);
    problem->M[n] = E - e * sin(E);
    problem->E_expect[n] = E;
  }
}

#endif
