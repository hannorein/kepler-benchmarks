#include <stdio.h>

#include "kepler_benchmarks/problems/basic.h"
#include "kepler_benchmarks/solvers/simple.h"

int main() {
  kb_problem* problem = kb_problem_alloc(100000);

  printf("ecc,runtime,size,error_max,error_mean\n");
  for (double e = 0.; e < 1; e += 0.001) {
    problem->e = e;
    kb_problem_setup_basic(problem);
    kb_summary result =
        kb_problem_run(problem, kb_solver_simple, kb_solver_simple_alloc, kb_solver_simple_free);
    printf("%.9f,", e);
    printf("%.9f,", result.runtime);
    printf("%ld,", result.size);
    printf("%.9f,", result.error_max);
    printf("%.9f\n", result.error_mean);
  }

  kb_problem_free(problem);

  return 0;
}
