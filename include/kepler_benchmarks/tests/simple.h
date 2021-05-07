#ifndef _KEPLER_BENCHMARKS_TESTS_SIMPLE_H_
#define _KEPLER_BENCHMARKS_TESTS_SIMPLE_H_

#include "test.h"

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
