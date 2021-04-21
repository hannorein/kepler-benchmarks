#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

struct kb_test {
    long N;
    double* M;
    double* e;
    double* E;
    double* E_solution;
};

struct kb_test_result {
    int pass;
    double runtime;
    double error_max;
    double error_mean;
};

struct kb_test kb_test_simple_new(long N){
    // Allocate memory
    struct kb_test test;
    test.N = N;
    test.M = malloc(sizeof(double)*N);
    test.e = malloc(sizeof(double)*N);
    test.E = malloc(sizeof(double)*N);
    test.E_solution = malloc(sizeof(double)*N);
    // Initialize test
    for (long i=0; i<N; i++){
        test.M[i] = i;
        test.e[i] = i;
        test.E_solution[i] = i;
    }
    return test;
}

void kb_test_simple_free(struct kb_test test){
    free(test.M);
    free(test.e);
    free(test.E);
    free(test.E_solution);
}

void kb_solver_simple(const double M, const double e, double* const E){
    *E = M*e;
}

struct kb_test_result kb_test_run(struct kb_test const test, void (*solver)(const double, const double, double* const)){
    struct timeval start, end;
    gettimeofday(&start, NULL);
    for (long i=0; i<test.N; i++){
        solver(test.M[i], test.e[i], &(test.E[i]));
    }
    gettimeofday(&end, NULL);
    struct kb_test_result result = {0};
    result.runtime = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
    for (long i=0; i<test.N; i++){
        double error = fabs(test.E[i]-test.E_solution[i]);
        result.error_max = max(result.error_max, error);
        result.error_mean += error/test.N;
    }
    return result;
}

int main(int argc, char* argv[]){
    long N = 1000000;
    struct kb_test test = kb_test_simple_new(N);
    struct kb_test_result result = kb_test_run(test, kb_solver_simple);
    printf("runtime:    %.9fs\n", result.runtime); 
    printf("error_max:  %.3e\n", result.error_max); 
    printf("error_mean: %.3e\n", result.error_mean); 
    kb_test_simple_free(test);
}
