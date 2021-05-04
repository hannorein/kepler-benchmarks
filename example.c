#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

//************************************
// Dummy Kepler solver
double kb_solver_simple(const double M, const double e, void* precalc){
    double precalculated_value = *(double*)precalc;
    // Return estimate for E
    return M * e * precalculated_value;
}

// Precalculate values
// Pointer can point to array, struct, etc
void* kb_solver_simple_alloc(const double e){
    double* p = malloc(sizeof(double));
    *p = sqrt(1.);
    return p;
}

// Free precalculated values
void kb_solver_simple_free(void* p){
    free(p);
}

struct kb_test_result {
    double runtime;
    double error_max;
    double error_mean;
};

//************************************
// Dummy Kepler test case
struct kb_test_result kb_test_simple(long N, double (*solver)(const double, const double, void* const), void* (*solver_alloc)(const double), void (*solver_free)(void* precalc)){
    // Allocate memory
    double* M = malloc(sizeof(double)*N);
    double* E = malloc(sizeof(double)*N);
    double* E_solution = malloc(sizeof(double)*N);
    // Initialize test
    for (long i=0; i<N; i++){
        M[i] = i;
        E_solution[i] = i;
    }
   
    struct timeval start, end;
    struct kb_test_result result = {0};
    // Test different eccentricities, might vary between tests
    for (double e=0.; e<1; e+=0.001){
        // Run solver setup
        void* precalc = NULL; // precalc can be zero
        if (solver_alloc){
            precalc = solver_alloc(e);
        }
        // Run test 
        gettimeofday(&start, NULL);
        for (long i=0; i<N; i++){
            E[i] = solver(M[i], e, precalc);
        }
        if (solver_free){
            solver_free(precalc);
        }
        gettimeofday(&end, NULL);
        result.runtime += (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
    }

    // Analyze results
    for (long i=0; i<N; i++){
        double error = fabs(E[i]-E_solution[i]);
        result.error_max = max(result.error_max, error);
        result.error_mean += error;
    }
    result.error_mean /= N;

    // Free memory
    free(M);
    free(E);
    free(E_solution);
    return result;
}

//////////////////////////////////////////////////////
int main(int argc, char* argv[]){
    struct kb_test_result result = kb_test_simple(100000, kb_solver_simple, kb_solver_simple_alloc, kb_solver_simple_free);
    // Solvers which do not precalculate anything can be called like this
    // struct kb_test_result result = kb_test_simple(100000, kb_solver_simple, NULL, NULL);
    printf("runtime:    %.9fs\n", result.runtime); 
    printf("error_max:  %.3e\n",  result.error_max); 
    printf("error_mean: %.3e\n",  result.error_mean); 
}
