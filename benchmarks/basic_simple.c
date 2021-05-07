#include <stdio.h>

#include "kepler_benchmarks/solvers/simple.h"
#include "kepler_benchmarks/problems/basic.h"

int main() {
  kb_test* test = kb_test_alloc(100000);

  printf("ecc,runtime,size,error_max,error_mean\n");
  for (double e = 0.; e < 1; e += 0.001) {
    kb_test_setup_simple(test, e);
    kb_test_result result =
        kb_test_run(test, kb_solver_simple, kb_solver_simple_alloc, kb_solver_simple_free);
    printf("%.9f,", e);
    printf("%.9f,", result.runtime);
    printf("%ld,", test->size);
    printf("%.9f,", result.error_max);
    printf("%.9f\n", result.error_mean);
  }

  kb_test_free(test);

  return 0;
}
