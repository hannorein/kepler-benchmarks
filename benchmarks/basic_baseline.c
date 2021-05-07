#include <stdio.h>

#include "kepler_benchmarks/problems/basic.h"
#include "kepler_benchmarks/solvers/baseline.h"

int main() {
  kb_problem* test = kb_problem_alloc(100000);

  printf("ecc,runtime,size,error_max,error_mean\n");
  for (double e = 0.; e < 1; e += 0.001) {
    kb_problem_setup_basic(test, e);
    kb_problem_result result = kb_problem_run(test, kb_solver_baseline, NULL, NULL);
    printf("%.9f,", e);
    printf("%.9f,", result.runtime);
    printf("%ld,", test->size);
    printf("%.9f,", result.error_max);
    printf("%.9f\n", result.error_mean);
  }

  kb_problem_free(test);

  return 0;
}
