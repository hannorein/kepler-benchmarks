#include "kepler_benchmarks/problems/basic.h"
#include "kepler_benchmarks/solvers.h"

int main() {
  kb_problem_run_ecc_grid(/*pre=*/0, /*size=*/100000, /*e_step=*/0.001,
                          /*setup=*/kb_problem_setup_basic,
                          /*solver=*/kb_solver_baseline, /*solver_alloc=*/NULL,
                          /*solver_free=*/NULL);
  return 0;
}
