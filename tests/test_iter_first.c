#include "helpers.h"
#include "kepler_benchmarks/solvers.h"

int main() { return run_solver_test(1e-8, kb_solver_iter_first, NULL, NULL); }
