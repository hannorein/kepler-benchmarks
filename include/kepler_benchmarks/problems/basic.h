#ifndef _KB_PROBLEMS_BASIC_H_
#define _KB_PROBLEMS_BASIC_H_

#include "kepler_benchmarks/problems/utils.h"

void kb_test_setup_simple(kb_test *test, double e) {
  test->fixed_ecc = 1;
  for (long n = 0; n < test->size; ++n) {
    double E = 2 * M_PI * n / (test->size);
    test->e[n] = e;
    test->M[n] = E - e * sin(E);
    test->E_expect[n] = E;
  }
}

#endif
